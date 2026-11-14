"""Modified nodal analysis, Newton-Raphson and transient integration."""

from __future__ import annotations

import copy
import math
import re
from dataclasses import dataclass, field, replace
from itertools import product

import numpy as np

from .parser import (
    CCCS,
    CCVS,
    GROUND_NAMES,
    VCCS,
    VCVS,
    AcAnalysis,
    BehavioralSource,
    Bjt,
    BjtModel,
    Capacitor,
    Circuit,
    DcAnalysis,
    Diode,
    DiodeModel,
    Inductor,
    ISource,
    MeasureSpec,
    Mosfet,
    MosModel,
    NoiseAnalysis,
    OpAnalysis,
    Resistor,
    StepAnalysis,
    TranAnalysis,
    VSource,
    _eval_expression,
    parse_value,
)
from .results import AcResult, DcSweep, NoiseResult, OperatingPoint, Transient

BOLTZMANN = 1.380649e-23
CHARGE = 1.602176634e-19


class ConvergenceError(RuntimeError):
    """Raised when Newton-Raphson fails to converge."""


@dataclass
class SolveOptions:
    reltol: float = 1e-3
    vntol: float = 1e-6
    abstol: float = 1e-9
    gmin: float = 1e-12
    maxiter: int = 100
    method: str = "trap"
    temp_c: float = 27.0


@dataclass
class _Context:
    mode: str
    time: float = 0.0
    dt: float = 0.0
    method: str = "trap"
    gmin: float = 1e-12
    source_scale: float = 1.0
    overrides: dict[str, float] = field(default_factory=dict)
    state: dict[str, dict[str, float]] = field(default_factory=dict)


class Simulation:
    """Run analyses for a parsed circuit."""

    def __init__(self, circuit: Circuit):
        self.circuit = circuit
        self.nodes = list(circuit.node_names)
        # Intrinsic BJT terminals, which may be hidden internal nodes when the
        # model declares series resistances (RB/RC/RE).
        self.internal_nodes: list[str] = []
        self.bjt_terminals: dict[str, tuple[str, str, str]] = {}
        self._bjt_parasitics: list[Resistor] = []
        self._build_bjt_parasitics()
        self.unknown_names = [*self.nodes, *self.internal_nodes]
        self.n_nodes = len(self.unknown_names)
        self.index = {name: i for i, name in enumerate(self.unknown_names)}
        self.branches = [
            el.name
            for el in circuit.elements
            if isinstance(el, (VSource, Inductor, VCVS, CCVS))
            or (isinstance(el, BehavioralSource) and el.kind == "v")
        ]
        self.branch_index = {
            name: self.n_nodes + i for i, name in enumerate(self.branches)
        }
        self.size = self.n_nodes + len(self.branches)
        # Devices whose current is computed from the solution (no branch row):
        # diodes and BJTs — exposed as I:<name> next to the branch currents.
        self.diodes = [el.name for el in circuit.elements if isinstance(el, Diode)]
        self.bjts = [el.name for el in circuit.elements if isinstance(el, Bjt)]
        self.behavioral = [
            el for el in circuit.elements if isinstance(el, BehavioralSource)
        ]
        self.opt = self._read_options()
        self._nonlinear = any(
            isinstance(el, (Diode, Mosfet, Bjt, BehavioralSource))
            for el in circuit.elements
        )
        self.mc_summary: dict[str, dict[str, float]] = {}
        self.yield_summary: dict[str, dict[str, float | None]] = {}
        self.yield_by_case: dict[str, dict[str, dict[str, float | None]]] = {}
        self.case_groups: list[dict] = []
        self.sensitivity_report: list[dict] = []

    # ------------------------------------------------------------- parasitics
    @staticmethod
    def _unique_node(base: str, used: set[str]) -> str:
        name = base
        suffix = 1
        while name in used or name in GROUND_NAMES:
            name = f"{base}_{suffix}"
            suffix += 1
        return name

    def _build_bjt_parasitics(self) -> None:
        """Create internal terminals and series resistors for BJT models."""
        used = set(self.nodes)
        for el in self.circuit.elements:
            if not isinstance(el, Bjt):
                continue
            model = self.circuit.models.get(el.model)
            if not isinstance(model, BjtModel):
                continue
            area = el.area if el.area > 0.0 else 1.0
            terminals: list[str] = []
            for tag, external, resistance in (
                ("c", el.nc, model.rc),
                ("b", el.nb, model.rb),
                ("e", el.ne, model.re),
            ):
                if resistance > 0.0:
                    name = self._unique_node(f"{el.name.lower()}#{tag}", used)
                    used.add(name)
                    self.internal_nodes.append(name)
                    self._bjt_parasitics.append(
                        Resistor(
                            f"{el.name}#r{tag}",
                            external,
                            name,
                            resistance / area,
                        )
                    )
                    terminals.append(name)
                else:
                    terminals.append(external)
            self.bjt_terminals[el.name] = (terminals[0], terminals[1], terminals[2])

    # --------------------------------------------------------------- options
    def _read_options(self) -> SolveOptions:
        raw = self.circuit.options

        def number(key: str, default: float) -> float:
            try:
                return float(parse_value(raw[key]))
            except (KeyError, ValueError):
                return default

        method = raw.get("method", "trap").lower()
        if method not in ("trap", "be"):
            method = "trap"
        return SolveOptions(
            reltol=number("reltol", 1e-3),
            vntol=number("vntol", 1e-6),
            abstol=number("abstol", 1e-9),
            gmin=number("gmin", 1e-12),
            maxiter=int(number("maxiter", 100)),
            method=method,
            temp_c=number("temp", 27.0),
        )

    def _thermal_voltage(self) -> float:
        return BOLTZMANN * (self.opt.temp_c + 273.15) / CHARGE

    def _solve_linear(self, A: np.ndarray, b: np.ndarray) -> np.ndarray:
        """Dense solve for small systems, SciPy sparse above a threshold."""
        if A.shape[0] >= 64:
            try:
                from scipy.sparse import csc_matrix
                from scipy.sparse.linalg import spsolve
            except ImportError:
                pass
            else:
                return np.asarray(spsolve(csc_matrix(A), b))
        return np.linalg.solve(A, b)

    # ------------------------------------------------------------- utilities
    def _node_voltage(self, nodes: dict[str, float], name: str) -> float:
        return 0.0 if name in GROUND_NAMES else nodes.get(name, 0.0)

    def _source_value(self, source: VSource | ISource, ctx: _Context) -> float:
        if source.name in ctx.overrides:
            value = ctx.overrides[source.name]
        elif ctx.mode == "tran":
            value = source.value_at(ctx.time)
        else:
            value = source.dc
        return value * ctx.source_scale

    def _find_source(self, name: str) -> VSource | ISource:
        for el in self.circuit.elements:
            if isinstance(el, (VSource, ISource)) and el.name.lower() == name.lower():
                return el
        raise ValueError(f"unknown source: {name}")

    def _branch_current_index(self, name: str) -> int:
        index = self.branch_index.get(name)
        if index is not None:
            return index
        for branch, idx in self.branch_index.items():
            if branch.lower() == name.lower():
                return idx
        raise ValueError(f"unknown controlling voltage source: {name}")

    # ------------------------------------------------------------ behavioral
    @staticmethod
    def _sanitize(name: str) -> str:
        return re.sub(r"\W", "_", name.strip().lower())

    def _runtime_params(self, x: np.ndarray, time: float = 0.0) -> dict[str, float]:
        params = dict(self.circuit.params)
        params["time"] = float(time)
        params["temper"] = float(self.opt.temp_c)
        params.setdefault("v__0", 0.0)
        params.setdefault("v__gnd", 0.0)
        for node, idx in self.index.items():
            params["v__" + self._sanitize(node)] = float(x[idx])
        for name, idx in self.branch_index.items():
            params["i__" + self._sanitize(name)] = float(x[idx])
        for name, current in self._diode_currents(x).items():
            params["i__" + self._sanitize(name)] = float(current)
        for name, current in self._bjt_currents(x).items():
            params["i__" + self._sanitize(name)] = float(current)
        return params

    def _behavioral_expr(self, el: BehavioralSource, x: np.ndarray) -> str:
        def voltage(match: re.Match[str]) -> str:
            first = "v__" + self._sanitize(match.group(1))
            if match.group(2):
                return f"({first} - v__{self._sanitize(match.group(2))})"
            return first

        expr = re.sub(
            r"\b[vV]\s*\(\s*([^,()]+?)\s*(?:,\s*([^()]+?)\s*)?\)",
            voltage,
            el.expr,
        )
        expr = re.sub(
            r"\b[iI]\s*\(\s*([^()]+?)\s*\)",
            lambda match: "i__" + self._sanitize(match.group(1)),
            expr,
        )
        return expr

    def _behavioral_value(
        self, el: BehavioralSource, x: np.ndarray, time: float = 0.0
    ) -> float:
        expr = self._behavioral_expr(el, x)
        return _eval_expression(
            expr, self._runtime_params(x, time), self.circuit.functions
        )

    def _behavioral_gradient(
        self, el: BehavioralSource, x: np.ndarray, time: float = 0.0
    ) -> np.ndarray:
        """Numerical d(expression)/d(state) for the AC linearization."""
        base = np.array(x, dtype=float)
        grad = np.zeros(self.size)
        for idx in range(self.size):
            step = 1e-6 * max(1.0, abs(base[idx]))
            plus = base.copy()
            minus = base.copy()
            plus[idx] += step
            minus[idx] -= step
            grad[idx] = (
                self._behavioral_value(el, plus, time)
                - self._behavioral_value(el, minus, time)
            ) / (2.0 * step)
        return grad

    # -------------------------------------------------------------- assembly
    def _assemble(self, x: np.ndarray, ctx: _Context) -> tuple[np.ndarray, np.ndarray]:
        n = self.size
        A = np.zeros((n, n))
        b = np.zeros(n)

        def stamp_g(i: int | None, j: int | None, g: float) -> None:
            if i is not None:
                A[i, i] += g
                if j is not None:
                    A[i, j] -= g
            if j is not None:
                A[j, j] += g
                if i is not None:
                    A[j, i] -= g

        def stamp_i(i: int | None, j: int | None, value: float) -> None:
            if i is not None:
                b[i] -= value
            if j is not None:
                b[j] += value

        def stamp_branch(k: int, i: int | None, j: int | None, value: float) -> None:
            if i is not None:
                A[i, k] += 1.0
            if j is not None:
                A[j, k] -= 1.0
            if i is not None:
                A[k, i] += 1.0
            if j is not None:
                A[k, j] -= 1.0
            b[k] += value

        for el in self.circuit.elements:
            if isinstance(el, Resistor):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                if el.value == 0:
                    raise ValueError(f"{el.name}: zero resistance")
                stamp_g(i, j, 1.0 / el.value)
            elif isinstance(el, Capacitor):
                if ctx.mode != "tran":
                    continue
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                prev = ctx.state.get(el.name, {"v": 0.0, "i": 0.0})
                if ctx.method == "be":
                    g = el.value / ctx.dt
                    ieq = -g * prev["v"]
                else:
                    g = 2.0 * el.value / ctx.dt
                    ieq = -g * prev["v"] - prev["i"]
                stamp_g(i, j, g)
                stamp_i(i, j, ieq)
            elif isinstance(el, Inductor):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                if ctx.mode == "tran":
                    prev = ctx.state.get(el.name, {"i": 0.0, "v": 0.0})
                    if ctx.method == "be":
                        req = el.value / ctx.dt
                        beq = -req * prev["i"]
                    else:
                        req = 2.0 * el.value / ctx.dt
                        beq = -req * prev["i"] - prev["v"]
                    if i is not None:
                        A[k, i] += 1.0
                    if j is not None:
                        A[k, j] -= 1.0
                    A[k, k] -= req
                    b[k] += beq
                else:
                    if i is not None:
                        A[k, i] += 1.0
                    if j is not None:
                        A[k, j] -= 1.0
                    b[k] += 0.0
                if i is not None:
                    A[i, k] += 1.0
                if j is not None:
                    A[j, k] -= 1.0
            elif isinstance(el, VSource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                stamp_branch(k, i, j, self._source_value(el, ctx))
            elif isinstance(el, ISource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                stamp_i(i, j, self._source_value(el, ctx))
            elif isinstance(el, BehavioralSource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                value = self._behavioral_value(el, x, ctx.time)
                grad = self._behavioral_gradient(el, x, ctx.time)
                if el.kind == "v":
                    k = self.branch_index[el.name]
                    stamp_branch(k, i, j, value - float(grad @ x))
                    A[k, :] -= grad
                else:
                    stamp_i(i, j, value)
                    shift = float(grad @ x)
                    if i is not None:
                        A[i, :] += grad
                        b[i] += shift
                    if j is not None:
                        A[j, :] -= grad
                        b[j] -= shift
            elif isinstance(el, VCVS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                stamp_branch(k, i, j, 0.0)
                nc1, nc2 = self.index.get(el.nc1), self.index.get(el.nc2)
                if nc1 is not None:
                    A[k, nc1] -= el.gain
                if nc2 is not None:
                    A[k, nc2] += el.gain
            elif isinstance(el, CCVS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                stamp_branch(k, i, j, 0.0)
                A[k, self._branch_current_index(el.control)] -= el.gain
            elif isinstance(el, VCCS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                nc1, nc2 = self.index.get(el.nc1), self.index.get(el.nc2)
                if i is not None:
                    if nc1 is not None:
                        A[i, nc1] += el.gain
                    if nc2 is not None:
                        A[i, nc2] -= el.gain
                if j is not None:
                    if nc1 is not None:
                        A[j, nc1] -= el.gain
                    if nc2 is not None:
                        A[j, nc2] += el.gain
            elif isinstance(el, CCCS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self._branch_current_index(el.control)
                if i is not None:
                    A[i, k] += el.gain
                if j is not None:
                    A[j, k] -= el.gain
            elif isinstance(el, Diode):
                self._stamp_diode(A, b, el, x)
            elif isinstance(el, Mosfet):
                self._stamp_mosfet(A, b, el, x)
            elif isinstance(el, Bjt):
                self._stamp_bjt(A, b, el, x, ctx)
            else:
                raise TypeError(f"unsupported element: {el!r}")

        for res in self._bjt_parasitics:
            i, j = self.index.get(res.n1), self.index.get(res.n2)
            stamp_g(i, j, 1.0 / res.value)

        for i in range(self.n_nodes):
            A[i, i] += ctx.gmin

        return A, b

    def _stamp_nonlinear(
        self,
        A: np.ndarray,
        b: np.ndarray,
        nodes: list[str],
        currents: list[float],
        jac: list[list[float]],
        x: np.ndarray,
    ) -> None:
        idx = [self.index.get(name) for name in nodes]
        vals = [0.0 if i is None else x[i] for i in idx]
        for row, ri in enumerate(idx):
            if ri is None:
                continue
            b[ri] += sum(jac[row][col] * vals[col] for col in range(len(nodes))) - currents[row]
            for col, ci in enumerate(idx):
                if ci is not None:
                    A[ri, ci] += jac[row][col]

    def _diode_saturation(self, model: DiodeModel) -> float:
        """Temperature-scaled saturation current of a diode model."""
        vt = self._thermal_voltage()
        temp_k = self.opt.temp_c + 273.15
        tnom_k = model.tnom + 273.15
        if temp_k == tnom_k or model.isat == 0.0:
            return model.isat
        ratio = temp_k / tnom_k
        return model.isat * ratio ** (model.xti / model.n) * math.exp(
            model.eg / (model.n * vt) * (ratio - 1.0)
        )

    def _diode_params(self, el: Diode, x: np.ndarray) -> tuple[float, float]:
        """Linearized (current, conductance) around the x operating point."""
        model = self.circuit.models.get(el.model)
        if not isinstance(model, DiodeModel):
            raise ValueError(f"{el.name}: unknown diode model {el.model!r}")
        i = self.index.get(el.n1)
        j = self.index.get(el.n2)
        v = (0.0 if i is None else x[i]) - (0.0 if j is None else x[j])
        vt = self._thermal_voltage()
        isat = self._diode_saturation(model)
        arg = min(v / (model.n * vt), 60.0)
        ex = math.exp(arg)
        current = el.area * isat * (ex - 1.0)
        conductance = el.area * isat / (model.n * vt) * ex
        return current, conductance

    def _stamp_diode(self, A: np.ndarray, b: np.ndarray, el: Diode, x: np.ndarray) -> None:
        current, conductance = self._diode_params(el, x)
        self._stamp_nonlinear(
            A,
            b,
            [el.n1, el.n2],
            [current, -current],
            [[conductance, -conductance], [-conductance, conductance]],
            x,
        )

    def _mosfet_state(self, el: Mosfet, x: np.ndarray) -> dict[str, float]:
        """Bias quantities of a MOSFET, normalized to its own polarity."""
        model = self.circuit.models.get(el.model)
        if not isinstance(model, MosModel):
            raise ValueError(f"{el.name}: unknown MOSFET model {el.model!r}")

        def voltage(name: str) -> float:
            idx = self.index.get(name)
            return 0.0 if idx is None else x[idx]

        sign = 1.0 if model.kind == "nmos" else -1.0
        vgs = sign * (voltage(el.ng) - voltage(el.ns))
        vds = sign * (voltage(el.nd) - voltage(el.ns))
        vbs = sign * (voltage(el.nb) - voltage(el.ns))

        vth = sign * model.vto
        dvth_dvbs = 0.0
        if model.gamma > 0.0:
            vbs_eff = min(vbs, 0.0)
            phib = max(model.phi - vbs_eff, 1e-3)
            vth += model.gamma * (math.sqrt(phib) - math.sqrt(model.phi))
            if vbs < 0.0:
                dvth_dvbs = -model.gamma / (2.0 * math.sqrt(phib))

        return {
            "sign": sign,
            "vgs": vgs,
            "vds": vds,
            "vbs": vbs,
            "vth": vth,
            "dvth_dvbs": dvth_dvbs,
            "beta": model.kp * (el.w / max(el.l, 1e-30)) * el.m,
            "vt": self._thermal_voltage(),
            "vov": vgs - vth,
        }

    def _mosfet_params(
        self, el: Mosfet, x: np.ndarray
    ) -> tuple[float, float, float, float]:
        """Signed (drain current, gm, gds, gmb) around the x operating point."""
        model = self.circuit.models.get(el.model)
        if not isinstance(model, MosModel):
            raise ValueError(f"{el.name}: unknown MOSFET model {el.model!r}")
        state = self._mosfet_state(el, x)
        sign = state["sign"]
        vds, vov, beta, vt = (
            state["vds"],
            state["vov"],
            state["beta"],
            state["vt"],
        )

        if vov <= 0.0:
            # Didactic subthreshold tail (NSUB on the model, default 1.5;
            # set NSUB=0 to disable it). Falls off exponentially with the
            # gate overdrive, slope nsub*Vt — ~89 mV/decade at 27C.
            if model.nsub > 0.0:
                vp = model.nsub * vt
                current_eff = 0.5 * beta * vp * vp * math.exp(vov / vp)
                gm = current_eff / vp
            else:
                current_eff = 0.0
                gm = 0.0
            gds = 0.0
        elif vds >= vov:
            f = 1.0 + model.lam * vds
            current_eff = 0.5 * beta * vov * vov * f
            gm = beta * vov * f
            gds = 0.5 * beta * vov * vov * model.lam
        else:
            f = 1.0 + model.lam * vds
            current_eff = beta * (vov * vds - 0.5 * vds * vds) * f
            gm = beta * vds * f
            gds = beta * (vov - vds) * f + beta * (
                vov * vds - 0.5 * vds * vds
            ) * model.lam

        gmb = -gm * state["dvth_dvbs"] if model.gamma > 0.0 else 0.0
        return sign * current_eff, gm, gds, gmb

    def _stamp_mosfet(self, A: np.ndarray, b: np.ndarray, el: Mosfet, x: np.ndarray) -> None:
        drain_current, gm, gds, gmb = self._mosfet_params(el, x)
        jd = [gds, gm, -(gm + gds + gmb), gmb]
        self._stamp_nonlinear(
            A,
            b,
            [el.nd, el.ng, el.ns, el.nb],
            [drain_current, 0.0, -drain_current, 0.0],
            [jd, [0.0, 0.0, 0.0, 0.0], [-v for v in jd], [0.0, 0.0, 0.0, 0.0]],
            x,
        )

    def _bjt_voltages(self, el: Bjt, x: np.ndarray) -> tuple[float, float, float]:
        """(polarity, Vbe, Vbc) normalized so an NPN is always positive."""
        model = self.circuit.models.get(el.model)
        if not isinstance(model, BjtModel):
            raise ValueError(f"{el.name}: unknown BJT model {el.model!r}")

        def voltage(name: str) -> float:
            idx = self.index.get(name)
            return 0.0 if idx is None else float(x[idx])

        nc, nb, ne = self.bjt_terminals.get(el.name, (el.nc, el.nb, el.ne))
        sign = 1.0 if model.kind == "npn" else -1.0
        return (
            sign,
            sign * (voltage(nb) - voltage(ne)),
            sign * (voltage(nb) - voltage(nc)),
        )

    def _bjt_charges(
        self, el: Bjt, x: np.ndarray
    ) -> tuple[float, BjtModel, float, float, float, float]:
        """Polarity, model, forward/reverse currents and their derivatives."""
        sign, vbe, vbc = self._bjt_voltages(el, x)
        model = self.circuit.models[el.model]
        assert isinstance(model, BjtModel)
        vt = self._thermal_voltage()
        isat = model.isat * el.area
        ef = math.exp(min(vbe / (model.nf * vt), 60.0))
        er = math.exp(min(vbc / (model.nr * vt), 60.0))
        forward = isat * (ef - 1.0)
        reverse = isat * (er - 1.0)
        df = isat / (model.nf * vt) * ef
        dr = isat / (model.nr * vt) * er
        return sign, model, forward, reverse, df, dr

    def _bjt_params(
        self, el: Bjt, x: np.ndarray
    ) -> tuple[float, float, float, list[list[float]]]:
        """(Ic, Ib, Ie, Jacobian) for an Ebers-Moll BJT."""
        model = self.circuit.models.get(el.model)
        if not isinstance(model, BjtModel):
            raise ValueError(f"{el.name}: unknown BJT model {el.model!r}")
        sign, vbe, vbc = self._bjt_voltages(el, x)
        vt = self._thermal_voltage()
        nfvt, nrvt = model.nf * vt, model.nr * vt
        ef = math.exp(min(vbe / nfvt, 60.0))
        er = math.exp(min(vbc / nrvt, 60.0))
        isat = model.isat * el.area
        bf = max(model.bf, 1e-12)
        br = max(model.br, 1e-12)
        forward = isat * (ef - 1.0)
        reverse = isat * (er - 1.0)
        df = isat / nfvt * ef
        dr = isat / nrvt * er
        transport = forward - reverse
        early = 1.0 - vbc / model.vaf if model.vaf > 0.0 else 1.0
        ic = transport * early - reverse / br

        # Gummel-Poon high-injection knee: extra base current at high Ic.
        ikf = model.ikf * el.area
        knee = 1.0 + forward / ikf if ikf > 0.0 else 1.0
        knee_slope = 1.0 + 2.0 * forward / ikf if ikf > 0.0 else 1.0
        ib = forward / bf * knee + reverse / br
        dib_dvbe = df / bf * knee_slope
        dib_dvbc = dr / br

        # Non-ideal (recombination) base currents.
        ise = model.ise * el.area
        isc = model.isc * el.area
        if ise > 0.0 and model.ne > 0.0:
            nevt = model.ne * vt
            ebe = math.exp(min(vbe / nevt, 60.0))
            ib += ise * (ebe - 1.0)
            dib_dvbe += ise / nevt * ebe
        if isc > 0.0 and model.nc > 0.0:
            ncvt = model.nc * vt
            ebc = math.exp(min(vbc / ncvt, 60.0))
            ib += isc * (ebc - 1.0)
            dib_dvbc += isc / ncvt * ebc

        ie = -(ic + ib)

        dic_dvbe = df * early
        dic_dvbc = (
            -dr * early
            - (transport / model.vaf if model.vaf > 0.0 else 0.0)
            - dr / br
        )
        # Rows are [Ic, Ib, Ie]; columns become [Vc, Vb, Ve].
        normalized = [
            [dic_dvbe, dic_dvbc],
            [dib_dvbe, dib_dvbc],
            [-(dic_dvbe + dib_dvbe), -(dic_dvbc + dib_dvbc)],
        ]
        jac = [
            [-row[1], row[0] + row[1], -row[0]]
            for row in normalized
        ]
        return sign * ic, sign * ib, sign * ie, jac

    def _stamp_bjt(
        self, A: np.ndarray, b: np.ndarray, el: Bjt, x: np.ndarray, ctx: _Context
    ) -> None:
        nc, nb, ne = self.bjt_terminals.get(el.name, (el.nc, el.nb, el.ne))
        ic, ib, ie, jac = self._bjt_params(el, x)
        self._stamp_nonlinear(A, b, [nc, nb, ne], [ic, ib, ie], jac, x)

        if ctx.mode != "tran":
            return

        _, model, forward, reverse, df, dr = self._bjt_charges(el, x)
        sign = 1.0 if model.kind == "npn" else -1.0
        area = el.area if el.area > 0.0 else 1.0
        idx = [self.index.get(node) for node in (nc, nb, ne)]
        values = [0.0 if i is None else float(x[i]) for i in idx]
        vbe = values[1] - values[2]
        vbc = values[1] - values[0]
        state = ctx.state.get(el.name, {})
        factor = 1.0 / ctx.dt if ctx.method == "be" else 2.0 / ctx.dt

        # Diffusion charges (TF/TR) and constant junction capacitances
        # (CJE/CJC). Newton linearizes the charge current around the last
        # iterate: i = i(v_k) + g·(v - v_k), with g = factor·dq/dv.
        charges = {
            "q_tf": sign * model.tf * forward,
            "q_tr": sign * model.tr * reverse,
            "q_cje": model.cje * area * vbe,
            "q_cjc": model.cjc * area * vbc,
        }
        slopes = {
            "g_tf": factor * model.tf * df,
            "g_tr": factor * model.tr * dr,
            "g_cje": factor * model.cje * area,
            "g_cjc": factor * model.cjc * area,
        }
        offsets = {}
        for key, charge in charges.items():
            current_key = "i" + key[1:]
            current = factor * (
                charge - state.get(key, 0.0)
            ) - state.get(current_key, 0.0)
            voltage = vbe if key in ("q_tf", "q_cje") else vbc
            offsets[key] = current - slopes["g" + key[1:]] * voltage

        i_be = offsets["q_tf"] + offsets["q_cje"]
        g_be = slopes["g_tf"] + slopes["g_cje"]
        i_bc = offsets["q_tr"] + offsets["q_cjc"]
        g_bc = slopes["g_tr"] + slopes["g_cjc"]
        self._stamp_cap_branch(A, b, idx[1], idx[2], i_be, g_be)
        self._stamp_cap_branch(A, b, idx[1], idx[0], i_bc, g_bc)

    @staticmethod
    def _stamp_cap_branch(
        A: np.ndarray,
        b: np.ndarray,
        i: int | None,
        j: int | None,
        current: float,
        conductance: float,
    ) -> None:
        if i is not None:
            b[i] -= current
            A[i, i] += conductance
            if j is not None:
                A[i, j] -= conductance
        if j is not None:
            b[j] += current
            A[j, j] += conductance
            if i is not None:
                A[j, i] -= conductance

    # ----------------------------------------------------------------- newton
    def _newton(self, ctx: _Context, x0: np.ndarray | None = None) -> tuple[np.ndarray, int]:
        x = np.zeros(self.size) if x0 is None else np.array(x0, dtype=float)
        if x.shape[0] != self.size:
            x = np.zeros(self.size)
        n_nodes = self.n_nodes

        for iteration in range(1, self.opt.maxiter + 1):
            A, b = self._assemble(x, ctx)
            try:
                xn = self._solve_linear(A, b)
            except np.linalg.LinAlgError as exc:
                raise ConvergenceError("singular matrix") from exc
            if not np.all(np.isfinite(xn)):
                raise ConvergenceError("non-finite solution")

            dx = xn - x
            if n_nodes and self._nonlinear:
                max_dv = float(np.max(np.abs(dx[:n_nodes])))
                if max_dv > 0.5:
                    xn = x + dx * (0.5 / max_dv)
                    dx = xn - x

            tol = np.empty(self.size)
            tol[:n_nodes] = self.opt.vntol + self.opt.reltol * np.maximum(
                np.abs(x[:n_nodes]), np.abs(xn[:n_nodes])
            )
            tol[n_nodes:] = self.opt.abstol + self.opt.reltol * np.maximum(
                np.abs(x[n_nodes:]), np.abs(xn[n_nodes:])
            )
            x = xn
            if np.all(np.abs(dx) <= tol):
                return x, iteration

        raise ConvergenceError(f"did not converge in {self.opt.maxiter} iterations")

    def _solve_dc(self, ctx: _Context, x0: np.ndarray | None) -> tuple[np.ndarray, int]:
        try:
            return self._newton(ctx, x0)
        except ConvergenceError as exc:
            direct_error = str(exc)

        # 1) gmin stepping at the full source level.
        x = x0
        gmin = 1e-2
        gmin_ok = True
        while gmin > self.opt.gmin * 10.0:
            try:
                x, _ = self._newton(replace(ctx, gmin=gmin), x)
            except ConvergenceError:
                gmin_ok = False
                break
            gmin /= 10.0
        if gmin_ok:
            try:
                return self._newton(replace(ctx, gmin=self.opt.gmin), x)
            except ConvergenceError:
                pass

        # 2) source stepping at the requested gmin.
        try:
            return self._source_stepping(ctx, x)
        except ConvergenceError:
            pass

        # 3) combined strategy: solve with a soft gmin, ramp the sources and
        #    then tighten gmin back to the requested value.
        soft_gmin = max(self.opt.gmin, 1e-6)
        try:
            x, iterations = self._source_stepping(
                replace(ctx, gmin=soft_gmin), x
            )
            while soft_gmin > self.opt.gmin * 10.0:
                soft_gmin = max(soft_gmin / 10.0, self.opt.gmin)
                x, iterations = self._newton(
                    replace(ctx, gmin=soft_gmin), x
                )
            return x, iterations
        except ConvergenceError as exc:
            diagnostics = self._residual_diagnostics(x, ctx)
            raise ConvergenceError(
                f"{direct_error}; source/gmin stepping failed; {diagnostics}"
            ) from exc

    def _source_stepping(
        self, ctx: _Context, x0: np.ndarray | None
    ) -> tuple[np.ndarray, int]:
        """Ramp independent sources from 0 to 1 with adaptive step control."""
        x = x0
        alpha = 0.0
        step = 0.1
        iterations = 0
        while alpha < 1.0 - 1e-12:
            target = min(alpha + step, 1.0)
            try:
                x, iterations = self._newton(
                    replace(ctx, source_scale=target), x
                )
            except ConvergenceError:
                step *= 0.5
                if step < 1e-4:
                    raise ConvergenceError(
                        f"source stepping stalled at alpha={alpha:.3f}"
                    )
                continue
            alpha = target
            step = min(step * 1.5, 0.25)
        assert x is not None
        return x, iterations

    def _residual_diagnostics(
        self, x: np.ndarray | None, ctx: _Context
    ) -> str:
        if x is None:
            vector = np.zeros(self.size)
        else:
            vector = np.asarray(x, dtype=float)
        if vector.shape[0] != self.size:
            return "no residual information"
        try:
            A, b = self._assemble(vector, ctx)
            residual = A @ vector - b
        except (ValueError, TypeError, FloatingPointError):
            return "no residual information"
        names = [*self.unknown_names, *self.branches]
        order = np.argsort(np.abs(residual))[::-1][:5]
        parts = [
            f"{names[index]}={residual[index]:.3g}"
            for index in order
            if index < len(names) and np.isfinite(residual[index])
        ]
        return "largest residuals: " + ", ".join(parts)


    # -------------------------------------------------------------- analyses
    def operating_point(self) -> OperatingPoint:
        ctx = _Context(mode="op", gmin=self.opt.gmin)
        x, iterations = self._solve_dc(ctx, None)
        return self._operating_point(x, iterations)

    # ----------------------------------------------------------- sensitivity
    @staticmethod
    def _op_signal(result: OperatingPoint, signal: str) -> float:
        """Value of ``V(node)``, ``V(a,b)`` or ``I(element)`` at an OP."""
        text = signal.strip()
        lowered = text.lower()
        if lowered.startswith("v(") and lowered.endswith(")"):
            parts = [part.strip() for part in text[2:-1].split(",")]
            if len(parts) > 2:
                raise ValueError(f"invalid voltage signal: {signal!r}")
            values = [
                (0.0 if name.lower() in GROUND_NAMES else result.nodes.get(name.lower(), 0.0))
                for name in parts
            ]
            return float(values[0]) if len(values) == 1 else float(values[0] - values[1])
        if lowered.startswith("i(") and lowered.endswith(")"):
            key = text[2:-1].strip()
            mapping = {name.lower(): value for name, value in result.currents.items()}
            if key.lower() not in mapping:
                raise ValueError(f"unknown current in .sens: {signal}")
            return float(mapping[key.lower()])
        mapping = {name.lower(): value for name, value in result.nodes.items()}
        if lowered in mapping:
            return float(mapping[lowered])
        currents = {name.lower(): value for name, value in result.currents.items()}
        if lowered in currents:
            return float(currents[lowered])
        raise ValueError(f"unknown signal in .sens: {signal}")

    @staticmethod
    def _sensitivity_targets(circuit: Circuit) -> list:
        """Elements whose value can be perturbed for sensitivity/MC."""
        return [
            el
            for el in circuit.elements
            if isinstance(el, Resistor)
            or (isinstance(el, (VSource, ISource)) and el.tol is not None)
        ]

    @staticmethod
    def _set_target_value(element, value: float) -> None:
        if isinstance(element, Resistor):
            element.value = value
        elif isinstance(element, (VSource, ISource)):
            element.dc = value

    @staticmethod
    def _target_value(element) -> float:
        return element.value if isinstance(element, Resistor) else element.dc

    def _sensitivity_at(
        self, circuit: Circuit, output: str, element, factor: float
    ) -> float:
        varied = copy.deepcopy(circuit)
        target = next(el for el in varied.elements if el.name == element.name)
        self._set_target_value(target, self._target_value(element) * factor)
        return self._op_signal(Simulation(varied).operating_point(), output)

    def sensitivity(
        self,
        output: str,
        top: int = 8,
        devices: int = 6,
        tol: float = 0.05,
        circuit: Circuit | None = None,
    ) -> dict:
        """Rank DC operating-point sensitivities and find the worst corners.

        The normalized sensitivity is ``(dY/dP) * P / Y`` so values with very
        different units can be compared. Only resistors and sources that
        declare ``TOL=`` perturb the DC solution; capacitors and inductors are
        open/short at DC and are therefore skipped.
        """
        base = circuit if circuit is not None else self.circuit
        nominal = self._op_signal(Simulation(base).operating_point(), output)
        entries: list[dict] = []
        for element in self._sensitivity_targets(base):
            value = self._target_value(element)
            if value == 0.0:
                continue
            relative = element.tol if element.tol is not None else tol
            if relative <= 0.0:
                continue
            # A small relative step gives an accurate derivative; ``relative``
            # is reserved for the worst-case corner sweep.
            step = 1e-4
            plus = self._sensitivity_at(base, output, element, 1.0 + step)
            minus = self._sensitivity_at(base, output, element, 1.0 - step)
            derivative = (plus - minus) / (2.0 * step * value)
            normalized = (
                derivative * value / nominal if abs(nominal) > 1e-12 else derivative
            )
            entries.append(
                {
                    "name": element.name,
                    "value": float(value),
                    "derivative": float(derivative),
                    "sensitivity": float(normalized),
                    "tol": float(relative),
                }
            )
        entries.sort(key=lambda entry: abs(entry["sensitivity"]), reverse=True)

        selected = entries[: min(max(devices, 1), 10)]
        worst = self._worst_case(base, output, selected)
        return {
            "output": output,
            "nominal": float(nominal),
            "entries": entries[: max(top, 1)],
            "ranked": entries,
            "worst_case": worst,
        }

    def _worst_case(self, circuit: Circuit, output: str, entries: list[dict]) -> dict:
        if not entries:
            nominal = self._op_signal(Simulation(circuit).operating_point(), output)
            return {"min": {"value": float(nominal), "case": {}}, "max": {"value": float(nominal), "case": {}}}
        factors = [(1.0 - entry["tol"], 1.0 + entry["tol"]) for entry in entries]
        best: dict | None = None
        worst: dict | None = None
        for combo in product(*factors):
            varied = copy.deepcopy(circuit)
            case: dict[str, float] = {}
            for entry, factor in zip(entries, combo):
                target = next(el for el in varied.elements if el.name == entry["name"])
                value = self._target_value(target) * factor
                self._set_target_value(target, value)
                case[entry["name"]] = float(value)
            try:
                value = self._op_signal(Simulation(varied).operating_point(), output)
            except ConvergenceError:
                continue
            point = {"value": float(value), "case": case}
            if best is None or value < best["value"]:
                best = point
            if worst is None or value > worst["value"]:
                worst = point
        if best is None or worst is None:
            nominal = self._op_signal(Simulation(circuit).operating_point(), output)
            return {"min": {"value": float(nominal), "case": {}}, "max": {"value": float(nominal), "case": {}}}
        return {"min": best, "max": worst}

    def _apply_sensitivity(self) -> None:
        self.sensitivity_report = []
        if not self.circuit.sensitivity:
            return
        if self.circuit.mc_runs > 0:
            self.circuit.warnings.append(
                ".sens is ignored with .mc; run sensitivity on the nominal circuit"
            )
            return
        if any(isinstance(a, StepAnalysis) for a in self.circuit.analyses):
            self.circuit.warnings.append(
                ".sens uses the base circuit and ignores .step cases"
            )
        for spec in self.circuit.sensitivity:
            try:
                self.sensitivity_report.append(
                    self.sensitivity(
                        spec.output,
                        top=spec.top,
                        devices=spec.devices,
                        tol=spec.tol,
                    )
                )
            except (ValueError, ConvergenceError) as exc:
                self.circuit.warnings.append(f".sens {spec.output}: {exc}")


    def _operating_point(self, x: np.ndarray, iterations: int) -> OperatingPoint:
        currents = {
            name: float(x[self.branch_index[name]]) for name in self.branches
        }
        currents.update(self._diode_currents(x))
        currents.update(self._bjt_currents(x))
        nodes = {name: float(x[i]) for name, i in self.index.items() if name in self.nodes}
        return OperatingPoint(
            nodes,
            currents,
            iterations,
            x.copy(),
            self._device_details(x),
        )

    def _device_details(self, x: np.ndarray) -> dict[str, dict[str, float | str]]:
        """Per-device bias report used by ``.op`` (and the GUI)."""
        details: dict[str, dict[str, float | str]] = {}

        def voltage(name: str) -> float:
            idx = self.index.get(name)
            return 0.0 if idx is None else float(x[idx])

        for el in self.circuit.elements:
            if isinstance(el, Diode):
                current, conductance = self._diode_params(el, x)
                details[el.name] = {
                    "type": "diode",
                    "model": el.model,
                    "region": "forward" if current > 0 else "reverse",
                    "vd": voltage(el.n1) - voltage(el.n2),
                    "id": current,
                    "gd": conductance,
                }
            elif isinstance(el, Mosfet):
                model = self.circuit.models.get(el.model)
                if not isinstance(model, MosModel):
                    continue
                state = self._mosfet_state(el, x)
                current, gm, gds, gmb = self._mosfet_params(el, x)
                if state["vov"] <= 0.0:
                    region = "subthreshold" if model.nsub > 0.0 else "cutoff"
                elif state["vds"] >= state["vov"]:
                    region = "saturation"
                else:
                    region = "linear"
                details[el.name] = {
                    "type": model.kind,
                    "model": model.name,
                    "region": region,
                    "id": current,
                    "vgs": state["vgs"],
                    "vds": state["vds"],
                    "vbs": state["vbs"],
                    "vth": state["vth"],
                    "vov": state["vov"],
                    "gm": gm,
                    "gds": gds,
                    "gmb": gmb,
                }
            elif isinstance(el, Bjt):
                model = self.circuit.models.get(el.model)
                if not isinstance(model, BjtModel):
                    continue
                _, vbe, vbc = self._bjt_voltages(el, x)
                ic, ib, ie, _ = self._bjt_params(el, x)
                if vbe < 0.4:
                    region = "cutoff"
                elif vbc < 0.0:
                    region = "forward active"
                else:
                    region = "saturation"
                details[el.name] = {
                    "type": model.kind,
                    "model": model.name,
                    "region": region,
                    "ic": ic,
                    "ib": ib,
                    "ie": ie,
                    "vbe": vbe,
                    "vbc": vbc,
                    "vce": vbe - vbc,
                    "beta": ic / ib if ib != 0.0 else float("nan"),
                }
        return details

    def _diode_currents(self, x: np.ndarray) -> dict[str, float]:
        """Current of every diode, computed from the solution vector."""
        by_name = {el.name: el for el in self.circuit.elements if isinstance(el, Diode)}
        return {name: self._diode_params(el, x)[0] for name, el in by_name.items()}

    def _bjt_currents(self, x: np.ndarray) -> dict[str, float]:
        """Collector current of every BJT (positive into the collector)."""
        by_name = {el.name: el for el in self.circuit.elements if isinstance(el, Bjt)}
        return {name: self._bjt_params(el, x)[0] for name, el in by_name.items()}

    def dc_sweep(self, source: str, start: float, stop: float, step: float) -> DcSweep:
        if step == 0:
            raise ValueError("DC sweep step must be non-zero")
        src = self._find_source(source)
        count = int(round((stop - start) / step))
        values = start + step * np.arange(max(count, 0) + 1)
        if step > 0:
            values = values[values <= stop + abs(step) * 1e-9]
        else:
            values = values[values >= stop - abs(step) * 1e-9]
        if values.size == 0:
            values = np.array([start])

        x: np.ndarray | None = None
        node_series: dict[str, list[float]] = {name: [] for name in self.nodes}
        current_series: dict[str, list[float]] = {
            name: [] for name in [*self.branches, *self.diodes, *self.bjts]
        }
        iterations: list[int] = []
        for value in values:
            ctx = _Context(
                mode="dc", overrides={src.name: float(value)}, gmin=self.opt.gmin
            )
            x, it = self._solve_dc(ctx, x)
            iterations.append(it)
            for name in self.nodes:
                node_series[name].append(float(x[self.index[name]]))
            for name in self.branches:
                current_series[name].append(float(x[self.branch_index[name]]))
            for name, current in self._diode_currents(x).items():
                current_series[name].append(current)
            for name, current in self._bjt_currents(x).items():
                current_series[name].append(current)

        return DcSweep(
            source=src.name,
            sweep=np.asarray(values, dtype=float),
            nodes={name: np.asarray(vals) for name, vals in node_series.items()},
            currents={name: np.asarray(vals) for name, vals in current_series.items()},
            iterations=iterations,
        )

    def _state_from_op(self, op: OperatingPoint) -> dict[str, dict[str, float]]:
        state: dict[str, dict[str, float]] = {}

        def voltage(name: str) -> float:
            if name in GROUND_NAMES:
                return 0.0
            if op.vector is not None:
                idx = self.index.get(name)
                if idx is not None:
                    return float(op.vector[idx])
            return op.nodes.get(name, 0.0)

        for el in self.circuit.elements:
            if isinstance(el, Capacitor):
                v = voltage(el.n1) - voltage(el.n2)
                state[el.name] = {"v": v, "i": 0.0}
            elif isinstance(el, Inductor):
                state[el.name] = {
                    "i": op.currents[el.name],
                    "v": voltage(el.n1) - voltage(el.n2),
                }
            elif isinstance(el, Bjt) and op.vector is not None:
                sign, model, forward, reverse, _, _ = self._bjt_charges(
                    el, op.vector
                )
                nc, nb, ne = self.bjt_terminals.get(
                    el.name, (el.nc, el.nb, el.ne)
                )
                area = el.area if el.area > 0.0 else 1.0
                vbe = voltage(nb) - voltage(ne)
                vbc = voltage(nb) - voltage(nc)
                state[el.name] = {
                    "q_tf": sign * model.tf * forward,
                    "q_tr": sign * model.tr * reverse,
                    "q_cje": model.cje * area * vbe,
                    "q_cjc": model.cjc * area * vbc,
                    "i_tf": 0.0,
                    "i_tr": 0.0,
                    "i_cje": 0.0,
                    "i_cjc": 0.0,
                }
        return state

    def _update_state(
        self,
        x: np.ndarray,
        previous: dict[str, dict[str, float]],
        dt: float,
        method: str,
    ) -> dict[str, dict[str, float]]:
        state: dict[str, dict[str, float]] = {}
        for el in self.circuit.elements:
            if isinstance(el, Capacitor):
                i = self.index.get(el.n1)
                j = self.index.get(el.n2)
                v = (0.0 if i is None else x[i]) - (0.0 if j is None else x[j])
                prev = previous.get(el.name, {"v": 0.0, "i": 0.0})
                if method == "be":
                    current = el.value / dt * (v - prev["v"])
                else:
                    current = 2.0 * el.value / dt * (v - prev["v"]) - prev["i"]
                state[el.name] = {"v": v, "i": current}
            elif isinstance(el, Inductor):
                i = self.index.get(el.n1)
                j = self.index.get(el.n2)
                current = x[self.branch_index[el.name]]
                v = (0.0 if i is None else x[i]) - (0.0 if j is None else x[j])
                state[el.name] = {"i": float(current), "v": v}
            elif isinstance(el, Bjt):
                sign, model, forward, reverse, _, _ = self._bjt_charges(el, x)
                nc, nb, ne = self.bjt_terminals.get(
                    el.name, (el.nc, el.nb, el.ne)
                )
                area = el.area if el.area > 0.0 else 1.0
                i = self.index.get(nb)
                j = self.index.get(ne)
                k = self.index.get(nc)
                vb = 0.0 if i is None else x[i]
                ve = 0.0 if j is None else x[j]
                vc = 0.0 if k is None else x[k]
                vbe = vb - ve
                vbc = vb - vc
                prev = previous.get(el.name, {})
                factor = 1.0 / dt if method == "be" else 2.0 / dt
                charges = {
                    "q_tf": sign * model.tf * forward,
                    "q_tr": sign * model.tr * reverse,
                    "q_cje": model.cje * area * vbe,
                    "q_cjc": model.cjc * area * vbc,
                }
                new: dict[str, float] = {}
                for key, charge in charges.items():
                    current_key = "i" + key[1:]
                    new[key] = charge
                    new[current_key] = factor * (
                        charge - prev.get(key, 0.0)
                    ) - prev.get(current_key, 0.0)
                state[el.name] = new
        return state

    def transient(
        self,
        tstep: float,
        tstop: float,
        tstart: float = 0.0,
        method: str | None = None,
        tmax: float | None = None,
        uic: bool = False,
    ) -> Transient:
        if tstep <= 0 or tstop <= 0:
            raise ValueError("transient times must be positive")
        method = method or self.opt.method
        if method not in ("trap", "be"):
            raise ValueError(f"unknown integration method: {method}")
        max_step = tstep if not tmax else min(tstep, tmax)

        if uic:
            x = np.zeros(self.size)
            for node, value in self.circuit.initial_conditions.items():
                if node in self.index:
                    x[self.index[node]] = value
            state = {}
            for el in self.circuit.elements:
                if isinstance(el, Capacitor):
                    if el.ic is not None:
                        v = el.ic
                        if el.n2 in GROUND_NAMES and el.n1 in self.index:
                            if el.n1 not in self.circuit.initial_conditions:
                                x[self.index[el.n1]] = v
                        elif el.n1 in GROUND_NAMES and el.n2 in self.index:
                            if el.n2 not in self.circuit.initial_conditions:
                                x[self.index[el.n2]] = -v
                    else:
                        i = self.index.get(el.n1)
                        j = self.index.get(el.n2)
                        v = (0.0 if i is None else x[i]) - (0.0 if j is None else x[j])
                    state[el.name] = {"v": v, "i": 0.0}
                elif isinstance(el, Inductor):
                    current = 0.0 if el.ic is None else el.ic
                    x[self.branch_index[el.name]] = current
                    i = self.index.get(el.n1)
                    j = self.index.get(el.n2)
                    v = (0.0 if i is None else x[i]) - (0.0 if j is None else x[j])
                    state[el.name] = {"i": current, "v": v}
            initial_iterations = 0
        else:
            op = self.operating_point()
            state = self._state_from_op(op)
            x = op.vector.copy() if op.vector is not None else np.zeros(self.size)
            initial_iterations = op.iterations

        time = 0.0
        dt = min(max_step, tstop)
        times = [0.0]
        rows = [x.copy()]
        iterations = [initial_iterations]

        while time < tstop - 1e-15:
            dt_try = min(dt, tstop - time)
            ctx = _Context(
                mode="tran",
                time=time + dt_try,
                dt=dt_try,
                method=method,
                gmin=self.opt.gmin,
                state=state,
            )
            try:
                xn, it = self._newton(ctx, x)
            except ConvergenceError:
                dt *= 0.5
                if dt < max(tstep * 1e-8, 1e-18):
                    raise
                continue

            time += dt_try
            x = xn
            state = self._update_state(x, state, dt_try, method)
            times.append(time)
            rows.append(x.copy())
            iterations.append(it)
            dt = min(max_step, dt_try * 1.5)

        arr = np.asarray(rows)
        t = np.asarray(times)
        if tstart > 0.0:
            mask = t >= tstart - 1e-18
            t = t[mask]
            arr = arr[mask]
            iterations = [it for it, keep in zip(iterations, mask) if keep]

        diode_series: dict[str, np.ndarray] = {}
        if self.diodes:
            by_name = {
                el.name: el for el in self.circuit.elements if isinstance(el, Diode)
            }
            diode_series = {
                name: np.asarray(
                    [self._diode_params(by_name[name], row)[0] for row in arr]
                )
                for name in self.diodes
            }
        bjt_series: dict[str, np.ndarray] = {}
        if self.bjts:
            by_name = {
                el.name: el for el in self.circuit.elements if isinstance(el, Bjt)
            }
            bjt_series = {
                name: np.asarray(
                    [self._bjt_params(by_name[name], row)[0] for row in arr]
                )
                for name in self.bjts
            }

        return Transient(
            time=t,
            nodes={name: arr[:, self.index[name]] for name in self.nodes},
            currents={
                **{name: arr[:, self.branch_index[name]] for name in self.branches},
                **diode_series,
                **bjt_series,
            },
            iterations=iterations,
        )

    # ------------------------------------------------------- small-signal AC
    @staticmethod
    def _ac_freqs(kind: str, points: int, fstart: float, fstop: float) -> np.ndarray:
        if fstart <= 0 or fstop <= 0 or fstart > fstop:
            raise ValueError(".ac needs 0 < fstart <= fstop")
        steps = max(int(points), 1)
        if kind == "lin":
            return np.linspace(fstart, fstop, steps)
        if fstop == fstart:
            return np.array([fstart])
        if kind == "dec":
            count = int(round(steps * math.log10(fstop / fstart)))
            return np.geomspace(fstart, fstop, count + 1)
        if kind == "oct":
            count = int(round(steps * math.log2(fstop / fstart)))
            return np.geomspace(fstart, fstop, count + 1)
        raise ValueError(f"unknown .ac kind: {kind!r}")

    def _assemble_ac(self, f: float, op: OperatingPoint) -> tuple[np.ndarray, np.ndarray]:
        """Complex small-signal MNA matrices at frequency f, linearized
        around the DC operating point `op`."""
        n = self.size
        A = np.zeros((n, n), dtype=complex)
        b = np.zeros(n, dtype=complex)
        s = 2j * math.pi * f
        x = op.vector

        def stamp_g(i: int | None, j: int | None, g: complex) -> None:
            if i is not None:
                A[i, i] += g
                if j is not None:
                    A[i, j] -= g
            if j is not None:
                A[j, j] += g
                if i is not None:
                    A[j, i] -= g

        def phasor(source: VSource | ISource) -> complex:
            mag, phase = source.ac if source.ac is not None else (0.0, 0.0)
            return mag * np.exp(1j * np.deg2rad(phase))

        for el in self.circuit.elements:
            if isinstance(el, Resistor):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                if el.value == 0:
                    raise ValueError(f"{el.name}: zero resistance")
                stamp_g(i, j, 1.0 / el.value)
            elif isinstance(el, Capacitor):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                stamp_g(i, j, s * el.value)
            elif isinstance(el, Inductor):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                if i is not None:
                    A[k, i] += 1.0
                    A[i, k] += 1.0
                if j is not None:
                    A[k, j] -= 1.0
                    A[j, k] -= 1.0
                A[k, k] -= s * el.value
            elif isinstance(el, VSource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                if i is not None:
                    A[i, k] += 1.0
                    A[k, i] += 1.0
                if j is not None:
                    A[j, k] -= 1.0
                    A[k, j] -= 1.0
                b[k] = phasor(el)     # no ac spec -> short circuit (0 V)
            elif isinstance(el, ISource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                value = phasor(el)
                if i is not None:
                    b[i] -= value
                if j is not None:
                    b[j] += value
            elif isinstance(el, BehavioralSource):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                grad = self._behavioral_gradient(el, op.vector)
                if el.kind == "v":
                    k = self.branch_index[el.name]
                    if i is not None:
                        A[i, k] += 1.0
                        A[k, i] += 1.0
                    if j is not None:
                        A[j, k] -= 1.0
                        A[k, j] -= 1.0
                    A[k, :] -= grad
                else:
                    if i is not None:
                        A[i, :] += grad
                    if j is not None:
                        A[j, :] -= grad
            elif isinstance(el, VCVS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                if i is not None:
                    A[i, k] += 1.0
                    A[k, i] += 1.0
                if j is not None:
                    A[j, k] -= 1.0
                    A[k, j] -= 1.0
                nc1, nc2 = self.index.get(el.nc1), self.index.get(el.nc2)
                if nc1 is not None:
                    A[k, nc1] -= el.gain
                if nc2 is not None:
                    A[k, nc2] += el.gain
            elif isinstance(el, CCVS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self.branch_index[el.name]
                if i is not None:
                    A[i, k] += 1.0
                    A[k, i] += 1.0
                if j is not None:
                    A[j, k] -= 1.0
                    A[k, j] -= 1.0
                A[k, self._branch_current_index(el.control)] -= el.gain
            elif isinstance(el, VCCS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                nc1, nc2 = self.index.get(el.nc1), self.index.get(el.nc2)
                if i is not None:
                    if nc1 is not None:
                        A[i, nc1] += el.gain
                    if nc2 is not None:
                        A[i, nc2] -= el.gain
                if j is not None:
                    if nc1 is not None:
                        A[j, nc1] -= el.gain
                    if nc2 is not None:
                        A[j, nc2] += el.gain
            elif isinstance(el, CCCS):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                k = self._branch_current_index(el.control)
                if i is not None:
                    A[i, k] += el.gain
                if j is not None:
                    A[j, k] -= el.gain
            elif isinstance(el, Diode):
                i, j = self.index.get(el.n1), self.index.get(el.n2)
                _, conductance = self._diode_params(el, x)
                stamp_g(i, j, conductance)
            elif isinstance(el, Mosfet):
                _, gm, gds, gmb = self._mosfet_params(el, x)
                jd = [gds, gm, -(gm + gds + gmb), gmb]
                jac = [jd, [0.0] * 4, [-v for v in jd], [0.0] * 4]
                idx = [self.index.get(node) for node in (el.nd, el.ng, el.ns, el.nb)]
                for row, ri in enumerate(idx):
                    if ri is None:
                        continue
                    for col, ci in enumerate(idx):
                        if ci is not None:
                            A[ri, ci] += jac[row][col]
            elif isinstance(el, Bjt):
                _, _, _, jac = self._bjt_params(el, x)
                nc, nb, ne = self.bjt_terminals.get(
                    el.name, (el.nc, el.nb, el.ne)
                )
                idx = [self.index.get(node) for node in (nc, nb, ne)]
                for row, ri in enumerate(idx):
                    if ri is None:
                        continue
                    for col, ci in enumerate(idx):
                        if ci is not None:
                            A[ri, ci] += jac[row][col]
                # Small-signal junction and transit-time capacitances.
                _, model, _, _, df, dr = self._bjt_charges(el, x)
                area = el.area if el.area > 0.0 else 1.0
                cbe = model.cje * area + model.tf * df
                cbc = model.cjc * area + model.tr * dr
                stamp_g(idx[1], idx[2], s * cbe)
                stamp_g(idx[1], idx[0], s * cbc)
            else:
                raise TypeError(f"unsupported element: {el!r}")

        for res in self._bjt_parasitics:
            i, j = self.index.get(res.n1), self.index.get(res.n2)
            stamp_g(i, j, 1.0 / res.value)

        for v in range(self.n_nodes):
            A[v, v] += self.opt.gmin
        return A, b

    def ac_sweep(
        self, kind: str, points: int, fstart: float, fstop: float
    ) -> AcResult:
        """Small-signal AC sweep, linearized around the DC operating point."""
        freqs = self._ac_freqs(kind, points, fstart, fstop)
        ctx = _Context(mode="op", gmin=self.opt.gmin)
        x, iterations = self._solve_dc(ctx, None)
        op = self._operating_point(x, iterations)
        nodes: dict[str, np.ndarray] = {
            name: np.empty(len(freqs), dtype=complex) for name in self.nodes
        }
        currents: dict[str, np.ndarray] = {
            name: np.empty(len(freqs), dtype=complex) for name in self.branches
        }
        for i, f in enumerate(freqs):
            A, b = self._assemble_ac(float(f), op)
            xac = self._solve_linear(A, b)
            for name in self.nodes:
                nodes[name][i] = xac[self.index[name]]
            for name, k in self.branch_index.items():
                currents[name][i] = xac[k]
        return AcResult(freqs, nodes, currents)

    def noise_sweep(
        self,
        output: str,
        input_source: str,
        kind: str,
        points: int,
        fstart: float,
        fstop: float,
    ) -> NoiseResult:
        """Output-referred noise using the AC transfer (adjoint) method."""
        freqs = self._ac_freqs(kind, points, fstart, fstop)
        ctx = _Context(mode="op", gmin=self.opt.gmin)
        x, iterations = self._solve_dc(ctx, None)
        op = self._operating_point(x, iterations)

        signal = output.strip()
        lowered = signal.lower()
        if lowered.startswith("v(") and lowered.endswith(")"):
            inner = signal[2:-1]
            parts = [part.strip() for part in inner.split(",")]
        else:
            parts = [signal]
        if len(parts) > 2:
            raise ValueError(f"invalid .noise output: {output!r}")
        output_vec = np.zeros(self.size)
        for node, factor in zip(parts, (1.0, -1.0)):
            idx = self.index.get(node.lower())
            if idx is not None:
                output_vec[idx] += factor

        # Unit excitation at the declared input source, used for input-referred
        # noise (gain from input to output).
        input_el = self._find_source(input_source)
        input_vec = np.zeros(self.size, dtype=complex)
        if isinstance(input_el, VSource):
            input_vec[self.branch_index[input_el.name]] = 1.0
        else:
            i, j = self.index.get(input_el.n1), self.index.get(input_el.n2)
            if i is not None:
                input_vec[i] -= 1.0
            if j is not None:
                input_vec[j] += 1.0

        kT = BOLTZMANN * (self.opt.temp_c + 273.15)
        sources: list[tuple[str, str, float]] = []
        for el in self.circuit.elements:
            if isinstance(el, Resistor) and el.value > 0:
                sources.append((el.n1, el.n2, 4.0 * kT / el.value))
            elif isinstance(el, Diode):
                current, _ = self._diode_params(el, x)
                sources.append((el.n1, el.n2, 2.0 * CHARGE * abs(current)))
            elif isinstance(el, Mosfet):
                _, gm, gds, _ = self._mosfet_params(el, x)
                sources.append((el.nd, el.ns, 4.0 * kT * (2.0 / 3.0) * abs(gm)))
                sources.append((el.nd, el.ns, 4.0 * kT * abs(gds)))
            elif isinstance(el, Bjt):
                ic, ib, _, _ = self._bjt_params(el, x)
                nc, nb, ne = self.bjt_terminals.get(
                    el.name, (el.nc, el.nb, el.ne)
                )
                sources.append((nc, ne, 2.0 * CHARGE * (abs(ic) + abs(ib))))
        for res in self._bjt_parasitics:
            if res.value > 0:
                sources.append((res.n1, res.n2, 4.0 * kT / res.value))

        psd = np.zeros(len(freqs), dtype=float)
        input_psd = np.full(len(freqs), np.nan, dtype=float)
        for i, f in enumerate(freqs):
            A, _ = self._assemble_ac(float(f), op)
            try:
                adjoint = np.linalg.solve(A.T, output_vec.astype(complex))
            except np.linalg.LinAlgError:
                adjoint = np.linalg.lstsq(A.T, output_vec.astype(complex), rcond=None)[0]
            for n1, n2, density in sources:
                i1, i2 = self.index.get(n1), self.index.get(n2)
                transfer = (0.0 if i1 is None else adjoint[i1]) - (
                    0.0 if i2 is None else adjoint[i2]
                )
                psd[i] += float(abs(transfer) ** 2) * density
            gain = complex(adjoint @ input_vec)
            if abs(gain) > 0.0:
                input_psd[i] = psd[i] / abs(gain) ** 2

        def integrate(values: np.ndarray) -> float:
            if len(freqs) < 2:
                return 0.0
            return float(
                np.nansum(0.5 * (values[1:] + values[:-1]) * np.diff(freqs))
            )

        integral = integrate(psd)
        input_integral = integrate(input_psd)
        return NoiseResult(
            freqs=freqs,
            output_psd=psd,
            output_name=parts[0],
            input_source=input_source,
            input_unit="V" if isinstance(input_el, VSource) else "A",
            total_rms=math.sqrt(max(integral, 0.0)),
            input_psd=input_psd,
            input_rms=math.sqrt(max(input_integral, 0.0)),
        )

    def _measure_series(
        self, result: DcSweep | Transient, signal: str
    ) -> tuple[np.ndarray, np.ndarray]:
        x = result.sweep if isinstance(result, DcSweep) else result.time
        text = signal.strip()
        lowered = text.lower()
        if lowered.startswith("v(") and lowered.endswith(")"):
            mapping, key = result.nodes, text[2:-1]
        elif lowered.startswith("i(") and lowered.endswith(")"):
            mapping, key = result.currents, text[2:-1]
        else:
            mapping, key = result.nodes, text
            if key.lower() not in {name.lower() for name in mapping}:
                mapping = result.currents
        lookup = {name.lower(): name for name in mapping}
        if key.lower() not in lookup:
            raise ValueError(f"unknown signal in .meas: {signal}")
        return x, mapping[lookup[key.lower()]]

    def _crossing(
        self,
        x: np.ndarray,
        y: np.ndarray,
        value: float,
        edge: str,
        occurrence: int,
    ) -> float:
        hits: list[float] = []
        for i in range(1, len(y)):
            y0, y1 = float(y[i - 1]), float(y[i])
            if edge == "rise":
                crossed = y0 < value <= y1
            else:
                crossed = y0 > value >= y1
            if crossed:
                fraction = (value - y0) / (y1 - y0) if y1 != y0 else 0.0
                hits.append(float(x[i - 1]) + fraction * (float(x[i]) - float(x[i - 1])))
        if len(hits) < occurrence:
            raise ValueError(
                f"signal never reaches {value:g} on a {edge} edge "
                f"(occurrence {occurrence})"
            )
        return hits[occurrence - 1]

    def _evaluate_measure(
        self, result: DcSweep | Transient, spec: MeasureSpec
    ) -> float:
        if spec.kind == "delay":
            assert spec.trig_signal is not None and spec.targ_signal is not None
            xt, yt = self._measure_series(result, spec.trig_signal)
            t0 = self._crossing(xt, yt, spec.trig_value, spec.trig_edge, spec.trig_occurrence)
            xr, yr = self._measure_series(result, spec.targ_signal)
            t1 = self._crossing(xr, yr, spec.targ_value, spec.targ_edge, spec.targ_occurrence)
            return t1 - t0

        assert spec.signal is not None
        x, y = self._measure_series(result, spec.signal)
        mask = np.ones_like(x, dtype=bool)
        if spec.from_x is not None:
            mask &= x >= spec.from_x
        if spec.to_x is not None:
            mask &= x <= spec.to_x
        x, y = x[mask], y[mask]
        if len(x) < 2:
            raise ValueError("measurement needs at least two points")

        if spec.kind in {"rise", "fall"}:
            span = float(np.max(y) - np.min(y))
            if span <= 0.0:
                raise ValueError("signal has no transition")
            low = spec.low if spec.low is not None else float(np.min(y)) + 0.1 * span
            high = spec.high if spec.high is not None else float(np.min(y)) + 0.9 * span
            if low >= high:
                raise ValueError("LOW must be below HIGH")
            edge = "rise" if spec.kind == "rise" else "fall"
            t_low = self._crossing(x, y, low, edge, spec.occurrence)
            t_high = self._crossing(x, y, high, edge, spec.occurrence)
            return (t_high - t_low) if spec.kind == "rise" else (t_low - t_high)
        if spec.kind == "cross":
            level = (
                spec.level
                if spec.level is not None
                else float(np.min(y) + np.max(y)) / 2.0
            )
            return self._crossing(x, y, level, spec.edge, spec.occurrence)
        if spec.kind == "width":
            level = (
                spec.level
                if spec.level is not None
                else float(np.min(y) + np.max(y)) / 2.0
            )
            start = self._crossing(x, y, level, spec.edge, spec.occurrence)
            other = "fall" if spec.edge == "rise" else "rise"
            after = x > start
            return (
                self._crossing(x[after], y[after], level, other, 1) - start
            )

        if spec.kind == "max":
            return float(np.max(y))
        if spec.kind == "min":
            return float(np.min(y))
        samples = y if spec.kind == "avg" else y * y
        value = float(np.sum(0.5 * (samples[1:] + samples[:-1]) * np.diff(x)))
        area = float(x[-1] - x[0])
        if area == 0.0:
            raise ValueError("measurement interval has zero length")
        if spec.kind == "rms":
            return float(math.sqrt(max(value / area, 0.0)))
        return float(value / area)

    def _ac_signal(self, result: AcResult, signal: str) -> np.ndarray:
        text = signal.strip()
        lowered = text.lower()
        if lowered.startswith("v(") and lowered.endswith(")"):
            mapping, key = result.nodes, text[2:-1]
        elif lowered.startswith("i(") and lowered.endswith(")"):
            mapping, key = result.currents, text[2:-1]
        else:
            mapping, key = result.nodes, text
            if key.lower() not in {name.lower() for name in mapping}:
                mapping = result.currents
        lookup = {name.lower(): name for name in mapping}
        if key.lower() not in lookup:
            raise ValueError(f"unknown signal in .meas: {signal}")
        return mapping[lookup[key.lower()]]

    def _evaluate_ac_measure(self, result: AcResult, spec: MeasureSpec) -> float:
        assert spec.signal is not None
        values = self._ac_signal(result, spec.signal)
        freqs = result.freqs
        mask = np.ones_like(freqs, dtype=bool)
        if spec.from_x is not None:
            mask &= freqs >= spec.from_x
        if spec.to_x is not None:
            mask &= freqs <= spec.to_x
        if not mask.any():
            raise ValueError("no AC point inside the measurement window")
        window = np.where(mask)[0]

        def index_at(value: float) -> int:
            return int(window[np.argmin(np.abs(freqs[window] - value))])

        if spec.kind == "gain":
            if spec.at is not None:
                index = index_at(spec.at)
            else:
                index = int(window[np.argmax(np.abs(values[window]))])
            return float(20.0 * np.log10(max(abs(values[index]), 1e-30)))
        if spec.kind == "phase":
            index = index_at(spec.at) if spec.at is not None else int(window[0])
            return float(np.rad2deg(np.angle(values[index])))
        if spec.kind == "loopgain":
            if spec.ref_signal is None:
                raise ValueError("LOOPGAIN needs an output and a reference signal")
            reference = self._ac_signal(result, spec.ref_signal)
            ratio = values / np.where(np.abs(reference) > 0.0, reference, 1e-30)
            if spec.at is not None:
                index = index_at(spec.at)
            else:
                index = int(window[np.argmax(np.abs(ratio[window]))])
            return float(20.0 * np.log10(max(abs(ratio[index]), 1e-30)))
        if spec.kind == "bw":
            gain = np.abs(values)
            low_gain = float(gain[window[0]])
            if low_gain <= 0.0:
                raise ValueError("zero AC gain")
            threshold = low_gain / math.sqrt(10.0 ** (spec.drop / 10.0))
            for position in range(1, len(window)):
                i0, i1 = window[position - 1], window[position]
                if gain[i1] <= threshold < gain[i0]:
                    g0, g1 = math.log10(gain[i0]), math.log10(gain[i1])
                    if g1 == g0:
                        return float(freqs[i1])
                    fraction = (math.log10(threshold) - g0) / (g1 - g0)
                    logf = math.log10(freqs[i0]) + fraction * (
                        math.log10(freqs[i1]) - math.log10(freqs[i0])
                    )
                    return float(10.0**logf)
            raise ValueError(f"gain never drops {spec.drop:g} dB")
        raise ValueError(f"unsupported AC measurement: {spec.kind}")

    def _apply_measurements(
        self, result: OperatingPoint | DcSweep | Transient | AcResult
    ) -> None:
        if isinstance(result, AcResult):
            kind = "ac"
        elif isinstance(result, DcSweep):
            kind = "dc"
        elif isinstance(result, Transient):
            kind = "tran"
        else:
            kind = "op"
        for spec in self.circuit.measures:
            if spec.analysis != kind:
                continue
            if isinstance(result, AcResult):
                try:
                    result.measurements[spec.name] = self._evaluate_ac_measure(
                        result, spec
                    )
                except ValueError as exc:
                    result.measurements[spec.name] = float("nan")
                    self.circuit.warnings.append(f".meas {spec.name}: {exc}")
                continue
            if not isinstance(result, (DcSweep, Transient)):
                result.measurements[spec.name] = float("nan")
                continue
            try:
                result.measurements[spec.name] = self._evaluate_measure(result, spec)
            except ValueError as exc:
                result.measurements[spec.name] = float("nan")
                self.circuit.warnings.append(f".meas {spec.name}: {exc}")

    def _apply_fourier(self, result: object) -> None:
        if not isinstance(result, Transient):
            return
        for spec in self.circuit.fourier:
            for signal in spec.signals:
                freqs, spectrum = result.fft(signal)
                dt = float(np.mean(np.diff(result.time)))
                resolution = 1.0 / (len(result.time) * dt)

                def amplitude(order: int) -> float:
                    target = order * spec.frequency
                    mask = np.abs(freqs - target) <= resolution
                    return float(np.max(spectrum[mask])) if mask.any() else 0.0

                label = re.sub(r"\W+", "_", signal).strip("_")
                for order in range(1, 6):
                    result.measurements[f"{label}_f{order}"] = amplitude(order)
                result.measurements[f"thd_{label}"] = result.thd(spec.frequency, signal)

    def _run_one(
        self, analysis: object
    ) -> OperatingPoint | DcSweep | Transient | AcResult | NoiseResult:
        if isinstance(analysis, OpAnalysis):
            result: object = self.operating_point()
        elif isinstance(analysis, DcAnalysis):
            result = self.dc_sweep(
                analysis.source, analysis.start, analysis.stop, analysis.step
            )
        elif isinstance(analysis, TranAnalysis):
            result = self.transient(
                analysis.tstep,
                analysis.tstop,
                analysis.tstart,
                self.opt.method,
                analysis.tmax,
                analysis.uic,
            )
        elif isinstance(analysis, AcAnalysis):
            result = self.ac_sweep(
                analysis.kind, analysis.points, analysis.fstart, analysis.fstop
            )
        elif isinstance(analysis, NoiseAnalysis):
            result = self.noise_sweep(
                analysis.output,
                analysis.input_source,
                analysis.kind,
                analysis.points,
                analysis.fstart,
                analysis.fstop,
            )
        else:
            raise TypeError(f"unsupported analysis: {analysis!r}")
        self._apply_measurements(result)
        self._apply_fourier(result)
        return result  # type: ignore[return-value]

    def _run_analyses(
        self, analyses: list, case: dict[str, float]
    ) -> list[OperatingPoint | DcSweep | Transient | AcResult | NoiseResult]:
        results = []
        for analysis in analyses:
            if isinstance(analysis, StepAnalysis):
                continue
            result = self._run_one(analysis)
            result.case = dict(case)
            results.append(result)
        return results

    def _summarize(self, results: list) -> dict[str, dict[str, float]]:
        samples: dict[str, list[float]] = {}
        for result in results:
            if isinstance(result, OperatingPoint):
                values = {
                    **{f"V({n})": v for n, v in result.nodes.items()},
                    **{f"I({n})": v for n, v in result.currents.items()},
                }
            elif isinstance(result, (DcSweep, Transient)):
                values = {
                    **{
                        f"V({n})": float(series[-1])
                        for n, series in result.nodes.items()
                    },
                    **{
                        f"I({n})": float(series[-1])
                        for n, series in result.currents.items()
                    },
                }
            elif isinstance(result, AcResult):
                values = {
                    f"|V({n})|": float(abs(series[-1]))
                    for n, series in result.nodes.items()
                }
            elif isinstance(result, NoiseResult):
                values = {"noise_out_rms": result.total_rms}
                if result.input_psd is not None:
                    values["noise_in_rms"] = result.input_rms
            else:
                values = {}
            for name, value in result.measurements.items():
                values[f"meas {name}"] = value
            for name, value in values.items():
                if np.isfinite(value):
                    samples.setdefault(name, []).append(float(value))
        summary: dict[str, dict[str, float]] = {}
        for name, values in samples.items():
            array = np.asarray(values)
            summary[name] = {
                "mean": float(np.mean(array)),
                "std": float(np.std(array, ddof=1)) if len(array) > 1 else 0.0,
                "min": float(np.min(array)),
                "max": float(np.max(array)),
            }
        return summary

    def _perturb(
        self,
        rng: np.random.Generator,
        global_tol: float,
        latent: float | None = None,
    ) -> None:
        """Apply one Monte Carlo perturbation draw to the circuit values.

        The default ``DIST=UNIFORM`` multiplies every device by
        ``1 + U(-tol, tol)``. ``DIST=GAUSS`` uses a Gaussian deviation with
        ``sigma = SIGMA`` (or ``tol / 3`` when ``SIGMA`` is omitted). ``RHO``
        mixes a common latent factor into every draw so nearby devices track
        each other while still varying individually.
        """
        gaussian = self.circuit.mc_dist == "gauss"
        rho = min(max(self.circuit.mc_rho, 0.0), 1.0)
        input_sigma = self.circuit.mc_sigma

        def deviation(tol: float) -> float:
            if not gaussian:
                return float(rng.uniform(-tol, tol))
            sigma = input_sigma if input_sigma is not None else tol / 3.0
            if rho > 0.0 and latent is not None:
                z = math.sqrt(rho) * latent + math.sqrt(1.0 - rho) * float(
                    rng.standard_normal()
                )
            else:
                z = float(rng.standard_normal())
            return sigma * z

        for el in self.circuit.elements:
            if isinstance(el, (Resistor, Capacitor, Inductor)):
                tol = global_tol if el.tol is None else el.tol
                if tol <= 0.0:
                    continue
                el.value *= max(1.0 + deviation(tol), 1e-6)
            elif isinstance(el, (VSource, ISource)) and el.tol:
                el.dc *= max(1.0 + deviation(el.tol), 1e-6)

    def _step_circuits(
        self, steps: list[StepAnalysis], base_analyses: list
    ) -> list[tuple[dict[str, float], Circuit, list]]:
        cases: list[tuple[dict[str, float], Circuit, list]] = []
        for values in product(*(step.values for step in steps)):
            case: dict[str, float] = {}
            overrides: dict[str, float] = {}
            temp: float | None = None
            source_values: dict[str, float] = {}
            for step, value in zip(steps, values):
                if step.kind == "param":
                    overrides[step.name] = value
                    case[step.name] = value
                elif step.kind == "temp":
                    temp = value
                    case["temp"] = value
                else:
                    source_values[step.name] = value
                    case[step.name] = value
            if overrides:
                from .parser import NetlistParser

                stepped = NetlistParser(
                    self.circuit.netlist_text,
                    self.circuit.base_dir or None,
                    overrides=overrides,
                    expr_seed=self.circuit.expr_seed,
                ).parse()
            else:
                stepped = copy.deepcopy(self.circuit)
            if temp is not None:
                stepped.options["temp"] = repr(temp)
            for name, value in source_values.items():
                for el in stepped.elements:
                    if (
                        isinstance(el, (VSource, ISource))
                        and el.name.lower() == name.lower()
                    ):
                        el.dc = value
                        break
            analyses = [
                a for a in stepped.analyses if not isinstance(a, StepAnalysis)
            ] or base_analyses
            cases.append((case, stepped, analyses))
        return cases

    def run(self) -> list[OperatingPoint | DcSweep | Transient | AcResult | NoiseResult]:
        steps = [a for a in self.circuit.analyses if isinstance(a, StepAnalysis)]
        self.mc_summary: dict[str, dict[str, float]] = {}
        base_analyses = [
            a for a in self.circuit.analyses if not isinstance(a, StepAnalysis)
        ]
        if steps:
            cases = self._step_circuits(steps, base_analyses)
        else:
            cases = [({}, self.circuit, base_analyses)]

        results: list = []
        if self.circuit.mc_runs > 0:
            global_tol = 0.05
            try:
                global_tol = float(parse_value(self.circuit.options.get("mc_tol", "0.05")))
            except ValueError:
                pass
            rng = np.random.default_rng(self.circuit.mc_seed)
            for case, circuit, analyses in cases:
                for run_index in range(1, self.circuit.mc_runs + 1):
                    varied = copy.deepcopy(circuit)
                    sim = Simulation(varied)
                    latent = (
                        float(rng.standard_normal())
                        if self.circuit.mc_rho > 0.0
                        else None
                    )
                    sim._perturb(rng, global_tol, latent)
                    results.extend(
                        sim._run_analyses(analyses, {**case, "mc": run_index})
                    )
            self.mc_summary = self._summarize(results)
            self._finalize_results(results)
            return results

        for case, circuit, analyses in cases:
            results.extend(Simulation(circuit)._run_analyses(analyses, case))
        self._finalize_results(results)
        return results

    @staticmethod
    def _case_label(case: dict[str, float]) -> str:
        parts = [f"{key}={value:g}" for key, value in case.items() if key != "mc"]
        return ", ".join(parts) if parts else "nominal"

    def _finalize_results(self, results: list) -> None:
        groups: dict[str, list[int]] = {}
        for index, result in enumerate(results):
            groups.setdefault(self._case_label(getattr(result, "case", {})), []).append(
                index
            )
        self.case_groups = [
            {"label": label, "indices": indices}
            for label, indices in groups.items()
        ]
        self._apply_sensitivity()
        if self.circuit.mc_runs <= 0:
            return
        self.yield_summary = self._measurement_stats(results)
        self.yield_by_case = {
            group["label"]: self._measurement_stats(
                [results[i] for i in group["indices"]]
            )
            for group in self.case_groups
        }

    def _measurement_stats(
        self, results: list
    ) -> dict[str, dict[str, float | None]]:
        """Mean/std/yield of every measured value across Monte Carlo runs."""
        limits = {
            spec.name: (spec.limit_min, spec.limit_max)
            for spec in self.circuit.measures
            if spec.limit_min is not None or spec.limit_max is not None
        }
        names = sorted(
            {
                name
                for result in results
                for name in getattr(result, "measurements", {})
            }
        )
        stats: dict[str, dict[str, float | None]] = {}
        for name in names:
            values = np.asarray(
                [result.measurements.get(name, np.nan) for result in results],
                dtype=float,
            )
            valid = values[np.isfinite(values)]
            low, high = limits.get(name, (None, None))
            if low is not None or high is not None:
                passed = valid[
                    (valid >= (low if low is not None else -np.inf))
                    & (valid <= (high if high is not None else np.inf))
                ]
            else:
                passed = None
            stats[name] = {
                "mean": float(np.mean(valid)) if valid.size else None,
                "std": float(np.std(valid)) if valid.size else None,
                "min": float(np.min(valid)) if valid.size else None,
                "max": float(np.max(valid)) if valid.size else None,
                "pass": int(passed.size) if passed is not None else None,
                "total": int(valid.size) if passed is not None else None,
                "yield": (
                    float(passed.size / valid.size)
                    if passed is not None and valid.size
                    else None
                ),
                "limit_min": low,
                "limit_max": high,
            }
        return stats
