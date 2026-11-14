"""Netlist parsing: values, devices, models, subcircuits and analyses."""

from __future__ import annotations

import ast
import math
import random
import re
from dataclasses import dataclass, field
from pathlib import Path

GROUND_NAMES = {"0", "gnd"}

_EXPR_CONSTANTS = {"pi": math.pi, "e": math.e, "true": 1.0, "false": 0.0}
_EXPR_FUNCTIONS = {
    "abs": abs,
    "acos": math.acos,
    "asin": math.asin,
    "atan": math.atan,
    "ceil": math.ceil,
    "cos": math.cos,
    "exp": math.exp,
    "floor": math.floor,
    "ln": math.log,
    "log": math.log,
    "log10": math.log10,
    "max": max,
    "min": min,
    "pow": pow,
    "sin": math.sin,
    "sqrt": math.sqrt,
    "tan": math.tan,
}
_BIN_OPS = {
    ast.Add: lambda a, b: a + b,
    ast.Sub: lambda a, b: a - b,
    ast.Mult: lambda a, b: a * b,
    ast.Div: lambda a, b: a / b,
    ast.FloorDiv: lambda a, b: a // b,
    ast.Mod: lambda a, b: a % b,
    ast.Pow: lambda a, b: a**b,
}
_UNARY_OPS = {
    ast.UAdd: lambda a: a,
    ast.USub: lambda a: -a,
    ast.Not: lambda a: 0.0 if a else 1.0,
}
_COMPARE_OPS = {
    ast.Eq: lambda a, b: a == b,
    ast.NotEq: lambda a, b: a != b,
    ast.Lt: lambda a, b: a < b,
    ast.LtE: lambda a, b: a <= b,
    ast.Gt: lambda a, b: a > b,
    ast.GtE: lambda a, b: a >= b,
}
# Numbers with a SPICE suffix inside an expression, e.g. "{2 * 1k}".
_EXPR_NUMBER = re.compile(
    r"(?<![\w.])((?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)([a-zA-Z]+)"
)

_MULTIPLIERS = {
    "t": 1e12,
    "g": 1e9,
    "k": 1e3,
    "m": 1e-3,
    "u": 1e-6,
    "n": 1e-9,
    "p": 1e-12,
    "f": 1e-15,
    "a": 1e-18,
}


class NetlistError(ValueError):
    """Raised when a netlist cannot be parsed."""


def _expand_expr_numbers(expr: str) -> str:
    """Turn SPICE suffixes inside expressions into plain Python numbers."""

    def repl(match: re.Match[str]) -> str:
        number = float(match.group(1))
        suffix = match.group(2).lower()
        if suffix.startswith("meg"):
            return repr(number * 1e6)
        if suffix.startswith("mil"):
            return repr(number * 25.4e-6)
        if suffix[0] in _MULTIPLIERS:
            return repr(number * _MULTIPLIERS[suffix[0]])
        return match.group(0)

    return _EXPR_NUMBER.sub(repl, expr)


def _eval_expression(
    expr: str,
    params: dict[str, float],
    functions: dict[str, tuple[list[str], str]] | None = None,
    rng: random.Random | None = None,
    depth: int = 0,
) -> float:
    """Evaluate the small arithmetic subset accepted inside ``{...}``.

    ``functions`` holds user ``.func`` definitions as ``(argument names,
    expression)`` pairs and ``rng`` enables the ``gauss(mean, sigma)``
    function. ``if(condition, then, else)`` is lazy: only the selected branch
    is evaluated.
    """
    if depth > 32:
        raise NetlistError("expression nesting too deep (recursive .func?)")
    # ``if`` is a Python keyword, so rename it to a callable before parsing.
    prepared = re.sub(r"\bif\s*\(", "__if__(", expr, flags=re.IGNORECASE)
    try:
        tree = ast.parse(_expand_expr_numbers(prepared), mode="eval")
    except SyntaxError as exc:
        raise NetlistError(f"invalid expression: {expr!r}") from exc

    def ev(node: ast.AST) -> float:
        if isinstance(node, ast.Expression):
            return ev(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            key = node.id.lower()
            if key in params:
                return float(params[key])
            if key in _EXPR_CONSTANTS:
                return float(_EXPR_CONSTANTS[key])
            raise NetlistError(f"unknown parameter in expression: {node.id}")
        if isinstance(node, ast.BinOp):
            op = _BIN_OPS.get(type(node.op))
            if op is None:
                raise NetlistError(f"unsupported operator in expression: {expr!r}")
            return float(op(ev(node.left), ev(node.right)))
        if isinstance(node, ast.UnaryOp):
            op = _UNARY_OPS.get(type(node.op))
            if op is None:
                raise NetlistError(f"unsupported operator in expression: {expr!r}")
            return float(op(ev(node.operand)))
        if isinstance(node, ast.Compare):
            left = ev(node.left)
            for operator, comparator in zip(node.ops, node.comparators):
                fn = _COMPARE_OPS.get(type(operator))
                if fn is None:
                    raise NetlistError(
                        f"unsupported comparison in expression: {expr!r}"
                    )
                right = ev(comparator)
                if not fn(left, right):
                    return 0.0
                left = right
            return 1.0
        if isinstance(node, ast.BoolOp):
            values = [ev(value) != 0.0 for value in node.values]
            if isinstance(node.op, ast.And):
                return 1.0 if all(values) else 0.0
            if isinstance(node.op, ast.Or):
                return 1.0 if any(values) else 0.0
            raise NetlistError(f"unsupported boolean operator: {expr!r}")
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            name = node.func.id.lower()
            if name in ("if", "__if__"):
                if len(node.args) != 3:
                    raise NetlistError("if() needs a condition, a then and an else")
                chosen = node.args[1] if ev(node.args[0]) != 0.0 else node.args[2]
                return ev(chosen)
            if name == "gauss":
                if rng is None:
                    raise NetlistError(
                        "gauss() is only available in .param/.func expressions"
                    )
                if len(node.args) > 2:
                    raise NetlistError("gauss() takes at most mean and sigma")
                mean = ev(node.args[0]) if node.args else 0.0
                sigma = ev(node.args[1]) if len(node.args) > 1 else 1.0
                return float(rng.gauss(mean, sigma))
            fn = _EXPR_FUNCTIONS.get(name)
            if fn is not None:
                return float(fn(*(ev(arg) for arg in node.args)))
            if functions and name in functions:
                arg_names, body = functions[name]
                if len(node.args) != len(arg_names):
                    raise NetlistError(
                        f".func {name} expects {len(arg_names)} argument(s)"
                    )
                local = dict(params)
                for arg_name, arg in zip(arg_names, node.args):
                    local[arg_name] = ev(arg)
                return _eval_expression(body, local, functions, rng, depth + 1)
            raise NetlistError(f"unsupported function in expression: {node.func.id}")
        raise NetlistError(f"unsupported expression: {expr!r}")

    return ev(tree)


def parse_value(
    text: str,
    params: dict[str, float] | None = None,
    functions: dict[str, tuple[list[str], str]] | None = None,
    rng: random.Random | None = None,
) -> float:
    """Parse a SPICE number such as ``10k``, ``2.2u`` or ``1e-3``.

    ``params`` enables ``{...}`` expressions and bare parameter names.
    """
    raw = text.strip()
    if not raw:
        raise NetlistError("empty value")
    if raw.startswith("{") and raw.endswith("}"):
        return _eval_expression(raw[1:-1], params or {}, functions, rng)
    s = raw.lower().replace("_", "")
    m = re.fullmatch(r"([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)([a-z]*)", s)
    if not m:
        if params and s in params:
            return float(params[s])
        raise NetlistError(f"invalid value: {text!r}")
    number = float(m.group(1))
    suffix = m.group(2)
    if suffix.startswith("meg"):
        return number * 1e6
    if suffix.startswith("mil"):
        return number * 25.4e-6
    if suffix and suffix[0] in _MULTIPLIERS:
        return number * _MULTIPLIERS[suffix[0]]
    return number


def _split_words(text: str) -> list[str]:
    """Split on whitespace/commas, keeping ``{...}`` and ``(...)`` intact."""
    parts: list[str] = []
    buffer: list[str] = []
    depth = 0
    for char in text.strip():
        if char in "({":
            depth += 1
        elif char in ")}":
            depth = max(0, depth - 1)
        if (char.isspace() or char == ",") and depth == 0:
            if buffer:
                parts.append("".join(buffer))
                buffer = []
        else:
            buffer.append(char)
    if buffer:
        parts.append("".join(buffer))
    return parts


def _args(text: str, params: dict[str, float] | None = None) -> list[float]:
    return [parse_value(part, params) for part in _split_words(text)]


@dataclass
class SinWave:
    offset: float
    amplitude: float
    freq: float
    delay: float = 0.0
    theta: float = 0.0
    phase: float = 0.0

    def __call__(self, t: float) -> float:
        if t < self.delay:
            return self.offset
        tt = t - self.delay
        return self.offset + self.amplitude * math.exp(-self.theta * tt) * math.sin(
            2.0 * math.pi * self.freq * tt + math.radians(self.phase)
        )


@dataclass
class PulseWave:
    v1: float
    v2: float
    td: float = 0.0
    tr: float = 0.0
    tf: float = 0.0
    pw: float = math.inf
    per: float = math.inf

    def __call__(self, t: float) -> float:
        if t < self.td:
            return self.v1
        tt = t - self.td
        if math.isfinite(self.per) and self.per > 0:
            tt = tt % self.per
        if tt < self.tr:
            return self.v1 + (self.v2 - self.v1) * (tt / self.tr if self.tr else 1.0)
        if tt < self.tr + self.pw:
            return self.v2
        if tt < self.tr + self.pw + self.tf:
            return self.v2 + (self.v1 - self.v2) * (
                (tt - self.tr - self.pw) / self.tf if self.tf else 1.0
            )
        return self.v1


@dataclass
class PwlWave:
    points: list[tuple[float, float]]

    def __call__(self, t: float) -> float:
        if t <= self.points[0][0]:
            return self.points[0][1]
        for (t0, v0), (t1, v1) in zip(self.points, self.points[1:]):
            if t <= t1:
                if t1 == t0:
                    return v1
                return v0 + (v1 - v0) * (t - t0) / (t1 - t0)
        return self.points[-1][1]


@dataclass
class Source:
    name: str
    n1: str
    n2: str
    dc: float = 0.0
    waveform: SinWave | PulseWave | PwlWave | None = None
    ac: tuple[float, float] | None = None
    tol: float | None = None

    def value_at(self, t: float) -> float:
        if self.waveform is not None:
            return float(self.waveform(t))
        return self.dc


@dataclass
class Resistor:
    name: str
    n1: str
    n2: str
    value: float
    tol: float | None = None


@dataclass
class Capacitor:
    name: str
    n1: str
    n2: str
    value: float
    ic: float | None = None
    tol: float | None = None


@dataclass
class Inductor:
    name: str
    n1: str
    n2: str
    value: float
    ic: float | None = None
    tol: float | None = None


@dataclass
class VSource(Source):
    pass


@dataclass
class ISource(Source):
    pass


@dataclass
class BehavioralSource:
    """Behavioral source ``B`` with a runtime expression.

    ``kind`` is ``"v"`` for an ideal voltage source and ``"i"`` for a
    current source. ``expr`` may reference ``V(node)``, ``I(element)`` and
    parameters from ``.param``.
    """

    name: str
    n1: str
    n2: str
    kind: str
    expr: str


@dataclass
class VCVS:
    name: str
    n1: str
    n2: str
    nc1: str
    nc2: str
    gain: float


@dataclass
class VCCS:
    name: str
    n1: str
    n2: str
    nc1: str
    nc2: str
    gain: float


@dataclass
class CCCS:
    name: str
    n1: str
    n2: str
    control: str
    gain: float


@dataclass
class CCVS:
    name: str
    n1: str
    n2: str
    control: str
    gain: float


@dataclass
class Diode:
    name: str
    n1: str
    n2: str
    model: str
    area: float = 1.0


@dataclass
class Mosfet:
    name: str
    nd: str
    ng: str
    ns: str
    nb: str
    model: str
    w: float = 1e-6
    l: float = 1e-6
    m: float = 1.0


@dataclass
class Bjt:
    name: str
    nc: str
    nb: str
    ne: str
    model: str
    area: float = 1.0


@dataclass
class DiodeModel:
    name: str
    isat: float = 1e-14
    n: float = 1.0
    rs: float = 0.0
    # Temperature scaling: IS(T) = IS(TNOM) * (T/TNOM)^(XTI/N)
    #                   * exp(EG/(N*Vt) * (T/TNOM - 1))
    xti: float = 3.0
    eg: float = 1.11
    tnom: float = 27.0


@dataclass
class MosModel:
    name: str
    kind: str = "nmos"
    vto: float = 0.0
    kp: float = 2e-5
    lam: float = 0.0
    gamma: float = 0.0
    phi: float = 0.6
    # Didactic subthreshold slope factor (NOT part of classic Level-1).
    # The exponential tail below the threshold decays with
    # nsub * Vt per e-fold: ~ nsub * Vt * ln(10) mV per decade (89 mV/dec
    # at nsub=1.5, T=27C). Set NSUB=0 in the model to disable the tail.
    nsub: float = 1.5


@dataclass
class BjtModel:
    name: str
    kind: str = "npn"
    isat: float = 1e-15
    bf: float = 100.0
    br: float = 1.0
    nf: float = 1.0
    nr: float = 1.0
    vaf: float = 0.0
    # Gummel-Poon subset: high-injection knee and non-ideal base currents.
    ikf: float = 0.0
    ise: float = 0.0
    ne: float = 1.5
    isc: float = 0.0
    nc: float = 2.0
    # Parasitic resistances and junction/transit-time capacitances.
    rb: float = 0.0
    rc: float = 0.0
    re: float = 0.0
    cje: float = 0.0
    cjc: float = 0.0
    tf: float = 0.0
    tr: float = 0.0


@dataclass
class OpAnalysis:
    pass


@dataclass
class DcAnalysis:
    source: str
    start: float
    stop: float
    step: float


@dataclass
class TranAnalysis:
    tstep: float
    tstop: float
    tstart: float = 0.0
    tmax: float | None = None
    uic: bool = False


@dataclass
class AcAnalysis:
    """Small-signal frequency sweep (.ac dec|oct|LIN points fstart fstop)."""

    kind: str  # "dec", "oct" or "lin"
    points: int
    fstart: float
    fstop: float


@dataclass
class StepAnalysis:
    """Parameter/temperature/source sweep (``.step``)."""

    kind: str          # "param", "temp" or "source"
    name: str          # parameter or source name (empty for temp)
    values: list[float]


@dataclass
class NoiseAnalysis:
    """Output/input-referred noise sweep (``.noise``)."""

    output: str        # e.g. "v(out)"
    input_source: str
    kind: str          # "dec", "oct" or "lin"
    points: int
    fstart: float
    fstop: float


@dataclass
class FourierSpec:
    """Harmonic analysis of a transient waveform (``.four``)."""

    frequency: float
    signals: list[str]


@dataclass
class MeasureSpec:
    """A small, useful subset of SPICE ``.meas`` statements."""

    analysis: str          # "op", "dc", "tran" or "ac"
    name: str
    kind: str              # "delay", "avg", "max", "min" or "rms"
    signal: str | None = None
    trig_signal: str | None = None
    trig_value: float = 0.0
    trig_edge: str = "rise"
    trig_occurrence: int = 1
    targ_signal: str | None = None
    targ_value: float = 0.0
    targ_edge: str = "rise"
    targ_occurrence: int = 1
    from_x: float | None = None
    to_x: float | None = None
    limit_min: float | None = None
    limit_max: float | None = None
    # Advanced transient/AC measurements.
    level: float | None = None
    low: float | None = None
    high: float | None = None
    at: float | None = None
    ref_signal: str | None = None
    drop: float = 3.0
    occurrence: int = 1
    edge: str = "rise"


@dataclass
class Subckt:
    name: str
    ports: list[str]
    lines: list[str]


@dataclass
class Circuit:
    title: str = ""
    elements: list = field(default_factory=list)
    models: dict[str, DiodeModel | MosModel | BjtModel] = field(default_factory=dict)
    analyses: list = field(default_factory=list)
    options: dict[str, str] = field(default_factory=dict)
    subckts: dict[str, Subckt] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    params: dict[str, float] = field(default_factory=dict)
    initial_conditions: dict[str, float] = field(default_factory=dict)
    measures: list[MeasureSpec] = field(default_factory=list)
    fourier: list[FourierSpec] = field(default_factory=list)
    functions: dict[str, tuple[list[str], str]] = field(default_factory=dict)
    mc_runs: int = 0
    mc_seed: int | None = None
    mc_dist: str = "uniform"
    mc_sigma: float | None = None
    mc_rho: float = 0.0
    expr_seed: int = 0
    netlist_text: str = ""
    base_dir: str = ""
    _node_order: dict[str, None] = field(default_factory=dict)

    def _add_node(self, node: str) -> None:
        if node not in GROUND_NAMES:
            self._node_order.setdefault(node, None)

    @property
    def node_names(self) -> list[str]:
        return list(self._node_order)


class NetlistParser:
    """Parse a SPICE-like netlist into a :class:`Circuit`."""

    def __init__(
        self,
        text: str,
        base_dir: str | Path | None = None,
        overrides: dict[str, float] | None = None,
        expr_seed: int = 0,
    ):
        self.text = text
        self.base_dir = Path(base_dir) if base_dir is not None else Path.cwd()
        self.circuit = Circuit()
        self._names: set[str] = set()
        self._main_lines: list[str] = []
        self._params: dict[str, float] = {}
        self._functions: dict[str, tuple[list[str], str]] = {}
        self._rng_seed = expr_seed
        self._rng = random.Random(expr_seed)
        self._overrides = {k.lower(): float(v) for k, v in (overrides or {}).items()}

    @classmethod
    def from_file(cls, path: str | Path) -> "NetlistParser":
        p = Path(path)
        return cls(p.read_text(encoding="utf-8"), base_dir=p.parent)

    def parse(self) -> Circuit:
        lines = self._preprocess(self.text, self.base_dir, 0)
        if lines and self._is_title(lines[0]):
            self.circuit.title = lines.pop(0).strip()
        self._collect_subckts(lines)
        all_lines = [
            *self._main_lines,
            *(line for sub in self.circuit.subckts.values() for line in sub.lines),
        ]
        self._scan_functions(all_lines)
        self._scan_params(all_lines)
        self.circuit.params.update(self._params)
        self.circuit.netlist_text = self.text
        self.circuit.base_dir = str(self.base_dir)
        self.circuit.expr_seed = self._rng_seed
        self._scan_models(self._main_lines)
        for sub in self.circuit.subckts.values():
            self._scan_models(sub.lines)
        for line in self._main_lines:
            self._parse_line(line)
        if not self.circuit.analyses:
            self.circuit.analyses.append(OpAnalysis())
        return self.circuit

    # ------------------------------------------------------------------ input
    def _preprocess(self, text: str, base_dir: Path, depth: int) -> list[str]:
        if depth > 16:
            raise NetlistError(".include nested too deeply")
        out: list[str] = []
        for raw in text.splitlines():
            line = raw.strip()
            if not line or line.startswith("*"):
                continue
            if ";" in line:
                line = line.split(";", 1)[0].strip()
            if not line:
                continue
            if line.startswith("+"):
                if not out:
                    raise NetlistError("continuation line without a previous line")
                out[-1] += " " + line[1:].strip()
                continue
            if line.lower().startswith(".include"):
                parts = line.split(None, 1)
                if len(parts) != 2:
                    raise NetlistError(".include needs a file path")
                path = (base_dir / parts[1].strip().strip("\"'")).expanduser()
                if not path.exists():
                    raise NetlistError(f"include not found: {path}")
                out.extend(self._preprocess(path.read_text(encoding="utf-8"), path.parent, depth + 1))
                continue
            out.append(line)
        return out

    @staticmethod
    def _is_title(line: str) -> bool:
        if line.startswith("."):
            return False
        toks = line.split()
        if len(toks) < 3:
            return True
        kind = toks[0][0].lower()
        if kind in "rclvi" and len(toks) >= 4:
            if kind in "vi" and toks[3].lower().startswith(
                ("sin(", "pulse(", "pwl(", "dc", "ac")
            ):
                return False
            if toks[3].startswith("{"):
                return False
            try:
                parse_value(toks[3])
                return False
            except NetlistError:
                pass
        if kind == "d" and len(toks) >= 4:
            return False
        if kind == "m" and len(toks) >= 5:
            return False
        if kind == "q" and len(toks) >= 5:
            return False
        if kind == "x" and len(toks) >= 4:
            return False
        if kind in "eg" and len(toks) >= 6:
            return False
        if kind in "fh" and len(toks) >= 5:
            return False
        if kind == "b" and len(toks) >= 4:
            return False
        return True

    def _collect_subckts(self, lines: list[str]) -> None:
        main: list[str] = []
        i = 0
        while i < len(lines):
            line = lines[i]
            low = line.lower()
            if low.startswith(".subckt"):
                parts = line.split()
                if len(parts) < 2:
                    raise NetlistError(".subckt needs a name")
                name = parts[1].lower()
                ports = [p.lower() for p in parts[2:]]
                body: list[str] = []
                i += 1
                while i < len(lines) and not lines[i].lower().startswith(".ends"):
                    body.append(lines[i])
                    i += 1
                if i >= len(lines):
                    raise NetlistError(f"unterminated .subckt {name}")
                self.circuit.subckts[name] = Subckt(name, ports, body)
            elif low.startswith(".ends"):
                raise NetlistError("unexpected .ends")
            else:
                main.append(line)
            i += 1
        self._main_lines = main

    def _scan_functions(self, lines: list[str]) -> None:
        """Collect ``.func name(args) = expression`` user functions."""
        for line in lines:
            if not line.lower().startswith(".func"):
                continue
            body = line[5:].strip()
            match = re.match(
                r"([A-Za-z_]\w*)\s*(?:\(([^)]*)\))?\s*=\s*(.+)", body
            )
            if not match:
                raise NetlistError(f"invalid .func: {line}")
            name = match.group(1).lower()
            args = [
                arg.strip().lower()
                for arg in (match.group(2) or "").split(",")
                if arg.strip()
            ]
            expr = match.group(3).strip()
            if expr.startswith("{") and expr.endswith("}"):
                expr = expr[1:-1]
            self._functions[name] = (args, expr)
        self.circuit.functions = dict(self._functions)

    def _scan_params(self, lines: list[str]) -> None:
        """Collect ``.param`` definitions, resolving references iteratively."""
        raw: dict[str, str] = {}
        for line in lines:
            if not line.lower().startswith(".param"):
                continue
            body = line[6:].strip()
            for match in re.finditer(
                r"([A-Za-z_]\w*)\s*=\s*(\{[^{}]*\}|[^\s,]+)", body
            ):
                raw[match.group(1).lower()] = match.group(2)

        pending = {name: text for name, text in raw.items() if name not in self._overrides}
        self._params.update(self._overrides)
        for _ in range(len(pending) + 2):
            changed = False
            for name, text in list(pending.items()):
                try:
                    value = parse_value(
                        text, self._params, self._functions, self._rng
                    )
                except NetlistError:
                    continue
                self._params[name] = value
                del pending[name]
                changed = True
            if not changed:
                break
        for name in pending:
            self.circuit.warnings.append(f"unresolved .param {name}")

    def _value(self, text: str) -> float:
        return parse_value(text, self._params, self._functions, self._rng)

    def _step_values(self, tokens: list[str]) -> list[float]:
        if not tokens:
            raise NetlistError(".step needs values")
        if tokens[0].lower() == "list":
            values = [self._value(token) for token in tokens[1:]]
            if not values:
                raise NetlistError(".step LIST needs at least one value")
            return values
        if len(tokens) < 3:
            raise NetlistError(".step needs start, stop and step")
        start, stop, step = (self._value(token) for token in tokens[:3])
        if step == 0:
            raise NetlistError(".step step must be non-zero")
        values: list[float] = []
        count = int(abs((stop - start) / step)) + 1
        for i in range(count + 1):
            value = start + i * step
            if step > 0 and value > stop + abs(step) * 1e-9:
                break
            if step < 0 and value < stop - abs(step) * 1e-9:
                break
            values.append(value)
        return values or [start]

    def _parse_step(self, line: str) -> StepAnalysis:
        tokens = _split_words(line)[1:]
        if not tokens:
            raise NetlistError(".step needs a target")
        first = tokens[0].lower()
        if first == "param":
            if len(tokens) < 4:
                raise NetlistError(".step param needs a name and values")
            return StepAnalysis("param", tokens[1].lower(), self._step_values(tokens[2:]))
        if first == "temp":
            return StepAnalysis("temp", "", self._step_values(tokens[1:]))
        return StepAnalysis("source", tokens[0], self._step_values(tokens[1:]))

    def _parse_noise(self, line: str) -> NoiseAnalysis:
        tokens = _split_words(line)[1:]
        if len(tokens) < 6:
            raise NetlistError(
                ".noise needs output, input source, dec/oct/lin, points, fstart and fstop"
            )
        kind = tokens[2].lower()
        if kind not in {"dec", "oct", "lin"}:
            raise NetlistError(f".noise kind must be dec, oct or lin (got {tokens[2]!r})")
        return NoiseAnalysis(
            output=tokens[0],
            input_source=tokens[1],
            kind=kind,
            points=int(float(self._value(tokens[3]))),
            fstart=self._value(tokens[4]),
            fstop=self._value(tokens[5]),
        )

    def _parse_initial_conditions(self, line: str) -> None:
        for match in re.finditer(
            r"[vV]\s*\(\s*([^)]+?)\s*\)\s*=\s*(\{[^{}]*\}|[^\s,]+)", line
        ):
            node = match.group(1).strip().lower()
            if node in GROUND_NAMES:
                continue
            self.circuit.initial_conditions[node] = self._value(match.group(2))
            self.circuit._add_node(node)

    def _parse_measure_condition(self, tokens: list[str]) -> dict:
        if not tokens:
            raise NetlistError(".meas condition needs a signal")
        result: dict = {"signal": tokens[0], "value": 0.0, "edge": "rise", "occurrence": 1}
        for token in tokens[1:]:
            if "=" not in token:
                raise NetlistError(f"invalid .meas condition: {' '.join(tokens)!r}")
            key, value = token.split("=", 1)
            key = key.lower()
            if key == "val":
                result["value"] = self._value(value)
            elif key == "rise":
                result["edge"] = "rise"
                result["occurrence"] = max(int(self._value(value)), 1)
            elif key == "fall":
                result["edge"] = "fall"
                result["occurrence"] = max(int(self._value(value)), 1)
            else:
                raise NetlistError(f"unsupported .meas keyword: {key}")
        return result

    def _parse_measure(self, line: str) -> None:
        tokens = _split_words(line)
        if len(tokens) < 4:
            raise NetlistError(".meas needs analysis, name and measurement type")
        analysis = tokens[1].lower()
        if analysis not in {"op", "dc", "tran", "ac"}:
            raise NetlistError(f"unsupported .meas analysis: {tokens[1]!r}")
        name = tokens[2]
        kind = tokens[3].upper()
        rest = tokens[4:]

        if kind == "TRIG":
            try:
                targ_at = next(
                    i for i, token in enumerate(rest) if token.upper() == "TARG"
                )
            except StopIteration as exc:
                raise NetlistError(".meas TRIG needs a TARG condition") from exc
            trig = self._parse_measure_condition(rest[:targ_at])
            targ = self._parse_measure_condition(rest[targ_at + 1 :])
            self.circuit.measures.append(
                MeasureSpec(
                    analysis=analysis,
                    name=name,
                    kind="delay",
                    trig_signal=trig["signal"],
                    trig_value=trig["value"],
                    trig_edge=trig["edge"],
                    trig_occurrence=trig["occurrence"],
                    targ_signal=targ["signal"],
                    targ_value=targ["value"],
                    targ_edge=targ["edge"],
                    targ_occurrence=targ["occurrence"],
                )
            )
            return

        transient_kinds = {"AVG", "MAX", "MIN", "RMS"}
        advanced_kinds = {"RISE", "FALL", "WIDTH", "CROSS"}
        ac_kinds = {"GAIN", "PHASE", "BW", "LOOPGAIN"}
        if kind not in transient_kinds | advanced_kinds | ac_kinds:
            raise NetlistError(f"unsupported .meas type: {tokens[3]!r}")
        if kind in ac_kinds and analysis != "ac":
            raise NetlistError(f".meas {kind} requires an .ac analysis")
        if kind in transient_kinds | advanced_kinds and analysis not in {"dc", "tran"}:
            raise NetlistError(f".meas {kind} requires a .dc or .tran analysis")
        if not rest:
            raise NetlistError(f".meas {kind} needs a signal")
        ref_signal = None
        if kind == "LOOPGAIN" and len(rest) > 1 and "=" not in rest[1]:
            ref_signal = rest[1]
            rest = [rest[0], *rest[2:]]
        from_x = to_x = None
        limit_min = limit_max = None
        level = low = high = at = None
        drop = 3.0
        occurrence = 1
        edge = "rise"
        for token in rest[1:]:
            if "=" not in token:
                raise NetlistError(f"invalid .meas argument: {token!r}")
            key, value = token.split("=", 1)
            key = key.lower()
            if key == "from":
                from_x = self._value(value)
            elif key == "to":
                to_x = self._value(value)
            elif key == "min":
                limit_min = self._value(value)
            elif key == "max":
                limit_max = self._value(value)
            elif key == "val":
                level = self._value(value)
            elif key == "low":
                low = self._value(value)
            elif key == "high":
                high = self._value(value)
            elif key == "at":
                at = self._value(value)
            elif key == "drop":
                drop = self._value(value)
            elif key == "rise":
                edge = "rise"
                occurrence = max(int(self._value(value)), 1)
            elif key == "fall":
                edge = "fall"
                occurrence = max(int(self._value(value)), 1)
            else:
                raise NetlistError(f"unsupported .meas keyword: {key}")
        if (
            limit_min is not None
            and limit_max is not None
            and limit_min > limit_max
        ):
            raise NetlistError(f".meas {name}: MIN is greater than MAX")
        self.circuit.measures.append(
            MeasureSpec(
                analysis=analysis,
                name=name,
                kind=kind.lower(),
                signal=rest[0],
                ref_signal=ref_signal,
                from_x=from_x,
                to_x=to_x,
                limit_min=limit_min,
                limit_max=limit_max,
                level=level,
                low=low,
                high=high,
                at=at,
                drop=drop,
                occurrence=occurrence,
                edge=edge,
            )
        )

    def _scan_models(self, lines: list[str]) -> None:
        for line in lines:
            if line.lower().startswith(".model"):
                self._parse_model(line)

    # ----------------------------------------------------------------- models
    def _parse_model(self, line: str) -> None:
        body = line[6:].strip()
        m = re.match(r"(\S+)\s+(.*)", body)
        if not m:
            raise NetlistError(f"invalid .model: {line}")
        name = m.group(1).lower()
        rest = m.group(2).replace("(", " ").replace(")", " ")
        toks = rest.split()
        if not toks:
            raise NetlistError(f"invalid .model: {line}")
        kind = toks[0].lower()
        params: dict[str, str] = {}
        for tok in toks[1:]:
            if "=" in tok:
                key, value = tok.split("=", 1)
                params[key.lower()] = value

        if kind == "d":
            self.circuit.models[name] = DiodeModel(
                name=name,
                isat=self._value(params.get("is", "1e-14")),
                n=self._value(params.get("n", "1")),
                rs=self._value(params.get("rs", "0")),
                xti=self._value(params.get("xti", "3")),
                eg=self._value(params.get("eg", "1.11")),
                tnom=self._value(params.get("tnom", "27")),
            )
        elif kind in ("nmos", "pmos"):
            self.circuit.models[name] = MosModel(
                name=name,
                kind=kind,
                vto=self._value(params.get("vto", "0")),
                kp=self._value(params.get("kp", "2e-5")),
                lam=self._value(params.get("lambda", params.get("lam", "0"))),
                gamma=self._value(params.get("gamma", "0")),
                phi=self._value(params.get("phi", "0.6")),
                nsub=self._value(params.get("nsub", "1.5")),
            )
        elif kind in ("npn", "pnp"):
            self.circuit.models[name] = BjtModel(
                name=name,
                kind=kind,
                isat=self._value(params.get("is", "1e-15")),
                bf=self._value(params.get("bf", "100")),
                br=self._value(params.get("br", "1")),
                nf=self._value(params.get("nf", "1")),
                nr=self._value(params.get("nr", "1")),
                vaf=self._value(params.get("vaf", "0")),
                ikf=self._value(params.get("ikf", "0")),
                ise=self._value(params.get("ise", "0")),
                ne=self._value(params.get("ne", "1.5")),
                isc=self._value(params.get("isc", "0")),
                nc=self._value(params.get("nc", "2")),
                rb=self._value(params.get("rb", "0")),
                rc=self._value(params.get("rc", "0")),
                re=self._value(params.get("re", "0")),
                cje=self._value(params.get("cje", "0")),
                cjc=self._value(params.get("cjc", "0")),
                tf=self._value(params.get("tf", "0")),
                tr=self._value(params.get("tr", "0")),
            )
        else:
            self.circuit.warnings.append(f"unsupported model type {kind!r}")

    # --------------------------------------------------------------- elements
    def _parse_line(self, line: str) -> None:
        low = line.lower()
        if low.startswith(".end"):
            return
        if low.startswith((".model", ".param", ".func")):
            return
        if low.startswith(".tran"):
            tokens = _split_words(line)[1:]
            uic = any(token.lower() == "uic" for token in tokens)
            parts = [token for token in tokens if token.lower() != "uic"]
            if len(parts) < 2:
                raise NetlistError(".tran needs tstep and tstop")
            tstart = self._value(parts[2]) if len(parts) > 2 else 0.0
            tmax = self._value(parts[3]) if len(parts) > 3 else None
            self.circuit.analyses.append(
                TranAnalysis(
                    self._value(parts[0]),
                    self._value(parts[1]),
                    tstart,
                    tmax,
                    uic,
                )
            )
        elif low.startswith(".dc"):
            parts = _split_words(line)[1:]
            if len(parts) < 4:
                raise NetlistError(".dc needs source, start, stop and step")
            if len(parts) > 4:
                self.circuit.warnings.append(
                    "nested .dc sweeps are not supported; using the first sweep"
                )
            self.circuit.analyses.append(
                DcAnalysis(
                    parts[0],
                    self._value(parts[1]),
                    self._value(parts[2]),
                    self._value(parts[3]),
                )
            )
        elif low.startswith(".ac"):
            parts = _split_words(line)[1:]
            if len(parts) < 4:
                raise NetlistError(".ac needs kind (dec/oct/lin), points, fstart and fstop")
            kind = parts[0].lower()
            if kind not in {"dec", "oct", "lin"}:
                raise NetlistError(f".ac kind must be dec, oct or lin (got {parts[0]!r})")
            if self._value(parts[2]) == 0:
                raise NetlistError(".ac fstart must be non-zero")
            self.circuit.analyses.append(
                AcAnalysis(
                    kind,
                    int(float(self._value(parts[1]))),
                    self._value(parts[2]),
                    self._value(parts[3]),
                )
            )
        elif low.startswith(".temp"):
            values = [self._value(token) for token in _split_words(line)[1:]]
            if not values:
                raise NetlistError(".temp needs a temperature")
            if len(values) > 1:
                self.circuit.warnings.append(
                    "multiple .temp values use the first; use '.step temp ...' "
                    "for a temperature sweep"
                )
            self.circuit.options["temp"] = repr(values[0])
        elif low.startswith(".step"):
            self.circuit.analyses.append(self._parse_step(line))
        elif low.startswith(".mc"):
            tokens = _split_words(line)[1:]
            if not tokens:
                raise NetlistError(".mc needs a number of runs")
            self.circuit.mc_runs = max(int(float(self._value(tokens[0]))), 0)
            for token in tokens[1:]:
                key, _, text = token.partition("=")
                key = key.lower()
                if key == "seed":
                    self.circuit.mc_seed = int(float(self._value(text)))
                elif key == "dist":
                    name = text.strip().lower()
                    if name in ("gauss", "gaussian", "normal"):
                        name = "gauss"
                    elif name in ("uniform", "flat"):
                        name = "uniform"
                    else:
                        raise NetlistError(
                            f"unsupported .mc distribution: {text!r}"
                        )
                    self.circuit.mc_dist = name
                elif key == "sigma":
                    sigma = self._value(text)
                    if sigma < 0.0:
                        raise NetlistError(".mc SIGMA must be non-negative")
                    self.circuit.mc_sigma = sigma
                elif key == "rho":
                    rho = self._value(text)
                    if not 0.0 <= rho <= 1.0:
                        raise NetlistError(".mc RHO must be between 0 and 1")
                    self.circuit.mc_rho = rho
                else:
                    raise NetlistError(f"unsupported .mc keyword: {key}")
            if self.circuit.mc_rho > 0.0 and self.circuit.mc_dist != "gauss":
                raise NetlistError(".mc RHO requires DIST=GAUSS")
        elif low.startswith(".noise"):
            self.circuit.analyses.append(self._parse_noise(line))
        elif low.startswith(".four"):
            tokens = _split_words(line)[1:]
            if len(tokens) < 2:
                raise NetlistError(".four needs a frequency and a signal")
            self.circuit.fourier.append(
                FourierSpec(self._value(tokens[0]), tokens[1:])
            )
        elif low.startswith(".ic"):
            self._parse_initial_conditions(line)
        elif low.startswith(".meas"):
            self._parse_measure(line)
        elif low.startswith(".options"):
            for tok in _split_words(line)[1:]:
                if "=" in tok:
                    key, value = tok.split("=", 1)
                    self.circuit.options[key.lower()] = value
        elif low == ".op" or low.startswith(".op "):
            self.circuit.analyses.append(OpAnalysis())
        elif low.startswith("."):
            if not low.startswith((".print", ".plot", ".save", ".width")):
                self.circuit.warnings.append(f"ignored directive: {line.split()[0]}")
        else:
            self._parse_element(line)

    def _parse_element(
        self,
        line: str,
        prefix: str | None = None,
        mapping: dict[str, str] | None = None,
        depth: int = 0,
    ) -> None:
        if depth > 32:
            raise NetlistError("subcircuit recursion too deep")
        low = line.lower().strip()
        if not low or low.startswith((".model", ".param", ".func", ".ends", ".subckt")):
            return
        toks = _split_words(line)
        if not toks:
            return
        base_name = toks[0]
        name = f"{prefix}:{base_name}" if prefix else base_name
        kind = base_name[0].lower()

        def node(token: str) -> str:
            n = token.lower()
            if n in GROUND_NAMES:
                return "0"
            if mapping and n in mapping:
                return mapping[n]
            if prefix:
                return f"{prefix}:{n}"
            return n

        def add(el, *nodes: str) -> None:
            if el.name in self._names:
                raise NetlistError(f"duplicate element name: {el.name}")
            self._names.add(el.name)
            for n in nodes:
                self.circuit._add_node(n)
            self.circuit.elements.append(el)

        if kind == "x":
            if len(toks) < 3:
                raise NetlistError(f"invalid subcircuit instance: {line}")
            sub_name = toks[-1].lower()
            sub = self.circuit.subckts.get(sub_name)
            if sub is None:
                raise NetlistError(f"unknown subcircuit: {sub_name}")
            nodes = [node(t) for t in toks[1:-1]]
            if len(nodes) != len(sub.ports):
                raise NetlistError(
                    f"{base_name}: expected {len(sub.ports)} nodes, got {len(nodes)}"
                )
            child_map = dict(zip(sub.ports, nodes))
            for subline in sub.lines:
                self._parse_element(subline, prefix=name, mapping=child_map, depth=depth + 1)
            return

        if kind in "rcl":
            if len(toks) < 4:
                raise NetlistError(f"invalid element: {line}")
            value = self._value(toks[3].split("=", 1)[-1])
            cls = {"r": Resistor, "c": Capacitor, "l": Inductor}[kind]
            n1, n2 = node(toks[1]), node(toks[2])
            initial = None
            tol = None
            for param in toks[4:]:
                key, _, text = param.partition("=")
                if key.lower() == "ic":
                    initial = self._value(text)
                elif key.lower() == "tol":
                    tol = self._value(text)
            if cls is Resistor:
                element = cls(name, n1, n2, value, tol)
            else:
                element = cls(name, n1, n2, value, initial, tol)
            add(element, n1, n2)
        elif kind in "vi":
            if len(toks) < 4:
                raise NetlistError(f"invalid source: {line}")
            n1, n2 = node(toks[1]), node(toks[2])
            tol = None
            source_tokens: list[str] = []
            for token in toks[3:]:
                if token.lower().startswith("tol="):
                    tol = self._value(token.split("=", 1)[1])
                else:
                    source_tokens.append(token)
            rest = " ".join(source_tokens)
            dc, waveform, ac = self._parse_source(rest)
            cls = VSource if kind == "v" else ISource
            add(cls(name, n1, n2, dc=dc, waveform=waveform, ac=ac, tol=tol), n1, n2)
        elif kind == "b":
            if len(toks) < 4:
                raise NetlistError(f"invalid behavioral source: {line}")
            n1, n2 = node(toks[1]), node(toks[2])
            match = re.match(r"\s*([vViI])\s*=\s*(\{.*\}|.+)$", " ".join(toks[3:]))
            if not match:
                raise NetlistError(f"invalid behavioral source: {line}")
            source_kind = match.group(1).lower()
            expr = match.group(2).strip()
            if expr.startswith("{") and expr.endswith("}"):
                expr = expr[1:-1]
            add(BehavioralSource(name, n1, n2, source_kind, expr), n1, n2)
        elif kind in "egfh":
            if kind in "eg":
                if len(toks) < 6:
                    raise NetlistError(f"invalid dependent source: {line}")
                n1, n2 = node(toks[1]), node(toks[2])
                nc1, nc2 = node(toks[3]), node(toks[4])
                gain = self._value(toks[5])
                cls = VCVS if kind == "e" else VCCS
                add(cls(name, n1, n2, nc1, nc2, gain), n1, n2, nc1, nc2)
            else:
                if len(toks) < 5:
                    raise NetlistError(f"invalid dependent source: {line}")
                n1, n2 = node(toks[1]), node(toks[2])
                control = f"{prefix}:{toks[3]}" if prefix else toks[3]
                gain = self._value(toks[4])
                cls = CCCS if kind == "f" else CCVS
                add(cls(name, n1, n2, control, gain), n1, n2)
        elif kind == "d":
            if len(toks) < 4:
                raise NetlistError(f"invalid diode: {line}")
            n1, n2 = node(toks[1]), node(toks[2])
            area = self._value(toks[4]) if len(toks) > 4 else 1.0
            add(Diode(name, n1, n2, toks[3].lower(), area), n1, n2)
        elif kind == "m":
            if len(toks) >= 6 and "=" not in toks[4] and "=" not in toks[5]:
                nd, ng, ns, nb = (node(t) for t in toks[1:5])
                model = toks[5].lower()
                params = toks[6:]
            elif len(toks) >= 5 and "=" not in toks[4]:
                nd, ng, ns = (node(t) for t in toks[1:4])
                nb = ns
                model = toks[4].lower()
                params = toks[5:]
            else:
                raise NetlistError(f"invalid MOSFET: {line}")
            w, l, mult = 1e-6, 1e-6, 1.0
            for param in params:
                if "=" not in param:
                    continue
                key, value = param.split("=", 1)
                key = key.lower()
                if key == "w":
                    w = self._value(value)
                elif key == "l":
                    l = self._value(value)
                elif key in ("m", "mult"):
                    mult = self._value(value)
            add(Mosfet(name, nd, ng, ns, nb, model, w, l, mult), nd, ng, ns, nb)
        elif kind == "q":
            if len(toks) < 5:
                raise NetlistError(f"invalid BJT: {line}")
            nc, nb, ne = (node(t) for t in toks[1:4])
            model = toks[4].lower()
            area = (
                self._value(toks[5])
                if len(toks) > 5 and "=" not in toks[5]
                else 1.0
            )
            add(Bjt(name, nc, nb, ne, model, area), nc, nb, ne)
        else:
            raise NetlistError(f"unsupported element: {line}")

    def _parse_source(
        self, text: str
    ) -> tuple[float, SinWave | PulseWave | PwlWave | None, tuple[float, float] | None]:
        dc = 0.0
        waveform: SinWave | PulseWave | PwlWave | None = None
        ac: tuple[float, float] | None = None

        m = re.match(r"([a-zA-Z]+)\s*\(([^)]*)\)(.*)", text, re.S)
        if m:
            func = m.group(1).lower()
            body = m.group(2)
            text = m.group(3).strip()
            if func == "sin":
                vals = _args(body, self._params)
                vals += [0.0] * (6 - len(vals))
                waveform = SinWave(*vals[:6])
            elif func == "pulse":
                vals = _args(body, self._params)
                if len(vals) < 2:
                    raise NetlistError("PULSE needs at least two values")
                vals += [0.0] * (7 - len(vals))
                if vals[6] <= 0:
                    vals[6] = math.inf
                if vals[5] <= 0:
                    vals[5] = math.inf
                waveform = PulseWave(*vals[:7])
            elif func == "pwl":
                vals = _args(body, self._params)
                if len(vals) < 2 or len(vals) % 2:
                    raise NetlistError("PWL needs time/value pairs")
                waveform = PwlWave(list(zip(vals[0::2], vals[1::2])))
            else:
                raise NetlistError(f"unsupported source function: {func}")

        toks = _split_words(text)
        i = 0
        if i < len(toks) and toks[i].lower() == "dc":
            i += 1
            if i >= len(toks):
                raise NetlistError("DC needs a value")
            dc = self._value(toks[i])
            i += 1
        elif i < len(toks):
            try:
                dc = self._value(toks[i])
                i += 1
            except NetlistError:
                pass
        if i < len(toks) and toks[i].lower() == "ac":
            mag = self._value(toks[i + 1]) if i + 1 < len(toks) else 1.0
            phase = self._value(toks[i + 2]) if i + 2 < len(toks) else 0.0
            ac = (mag, phase)
        return dc, waveform, ac
