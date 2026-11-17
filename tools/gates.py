"""Compound-gate registry, testbench composer and truth-table helpers.

Each cell lives in ``cells/sky130/<name>.cir`` as a SPICE ``.subckt`` with a
``vdd`` port and an output port ``y``.  The expected logic is stored here so
that both the pytest suite and the ``verify_gates.py`` script share one source
of truth.

Multi-output cells (e.g. ``fa`` with ``sum``/``cout``) register one ``Gate`` per
output: the entries share ``cell``/``inputs`` and differ in ``output``, so the
single-output test loop keeps working unchanged.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CELLS = ROOT / "cells" / "sky130"
MODELS = CELLS / "models.spice"

VDD = 1.8  # sky130 core supply
LOW, HIGH = 0.25 * VDD, 0.75 * VDD  # interpretation margins


@dataclass
class Gate:
    name: str
    inputs: tuple[str, ...]
    logic: str  # human-readable expression
    expect: callable  # inputs (0/1) -> 0/1
    output: str = "y"  # output node under test
    cell: str = ""  # netlist basename, defaults to ``name``

    def __post_init__(self) -> None:
        if not self.cell:
            self.cell = self.name

    @property
    def port_order(self) -> tuple[str, ...]:
        # X-instance port order must match the .subckt declaration.
        return dict_port_order[self.cell]

    @property
    def outputs(self) -> tuple[str, ...]:
        return dict_output_ports.get(self.cell, ("y",))


# Port order as declared in each netlist file (.subckt <name> ... y vdd).
dict_port_order: dict[str, tuple[str, ...]] = {
    "inv": ("a",),
    "nand2": ("a", "b"),
    "nor2": ("a", "b"),
    "aoi21": ("a", "b", "c"),
    "oai21": ("a", "b", "c"),
    "aoi22": ("a", "b", "c", "d"),
    "oai22": ("a", "b", "c", "d"),
    "mux21": ("i0", "i1", "s"),
    "fa": ("a", "b", "cin"),
}

# Output ports after the inputs (default: the single port `y`).
dict_output_ports: dict[str, tuple[str, ...]] = {
    "fa": ("sum", "cout"),
}

GATES: dict[str, Gate] = {
    "inv": Gate(
        "inv", ("a",), "y = !a", lambda a: 1 - a
    ),
    "nand2": Gate(
        "nand2", ("a", "b"), "y = !(a & b)", lambda a, b: 1 - (a & b)
    ),
    "nor2": Gate(
        "nor2", ("a", "b"), "y = !(a | b)", lambda a, b: 1 - (a | b)
    ),
    "aoi21": Gate(
        "aoi21", ("a", "b", "c"), "y = !((a & b) | c)",
        lambda a, b, c: 1 - ((a & b) | c),
    ),
    "oai21": Gate(
        "oai21", ("a", "b", "c"), "y = !((a | b) & c)",
        lambda a, b, c: 1 - ((a | b) & c),
    ),
    "aoi22": Gate(
        "aoi22", ("a", "b", "c", "d"), "y = !((a & b) | (c & d))",
        lambda a, b, c, d: 1 - ((a & b) | (c & d)),
    ),
    "oai22": Gate(
        "oai22", ("a", "b", "c", "d"), "y = !((a | b) & (c | d))",
        lambda a, b, c, d: 1 - ((a | b) & (c | d)),
    ),
    "mux21": Gate(
        "mux21", ("i0", "i1", "s"), "y = !(i0 if s else i1)",
        lambda i0, i1, s: 1 - (i0 if s else i1),
    ),
    # Full adder: one netlist (fa.cir), one registry row per output.
    "fa_sum": Gate(
        "fa_sum", ("a", "b", "cin"), "sum = a ^ b ^ cin",
        lambda a, b, cin: a ^ b ^ cin, output="sum", cell="fa",
    ),
    "fa_cout": Gate(
        "fa_cout", ("a", "b", "cin"), "cout = (a & b) | (cin & (a ^ b))",
        lambda a, b, cin: (a & b) | (cin & (a ^ b)), output="cout", cell="fa",
    ),
}


def combos(gate: Gate) -> list[tuple[int, ...]]:
    """All input combinations for a gate, MSB first."""
    return list(product((0, 1), repeat=len(gate.inputs)))


def decode(voltage: float) -> int:
    """Interpret an output voltage as a logic level."""
    if voltage <= LOW:
        return 0
    if voltage >= HIGH:
        return 1
    raise ValueError(f"ambiguous output voltage: {voltage:.4g} V")


def build_testbench(gate: Gate, inputs: tuple[int, ...]) -> str:
    """Compose a self-contained `.op` netlist that drives the cell."""
    lines = [
        f"* testbench: {gate.name} {inputs}",
        f".include {MODELS}",
        f".include {CELLS / (gate.cell + '.cir')}",
        f"VDD vdd 0 {VDD}",
    ]
    for name, level in zip(gate.inputs, inputs):
        lines.append(f"V{name} {name} 0 {level * VDD}")
    ports = " ".join(dict_port_order[gate.cell])
    outs = " ".join(gate.outputs)
    lines.append(f"XG {ports} {outs} vdd {gate.cell}")
    lines.append(".op")
    lines.append(".end")
    return "\n".join(lines) + "\n"


def to_docstring() -> str:
    return "\n".join(f"{g.name}: {g.logic}" for g in GATES.values())
