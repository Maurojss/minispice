#!/usr/bin/env python
"""Print and verify the truth table of every sky130 compound cell.

Usage::

    python scripts/verify_gates.py [gate ...]

Exits non-zero if any output does not match the expected logic.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

from minispice import NetlistParser, Simulation
from tools.gates import GATES, LOW, VDD, build_testbench, combos


def run_case(gate, inputs):
    circuit = NetlistParser(build_testbench(gate, inputs)).parse()
    y = Simulation(circuit).operating_point().nodes["y"]
    got = 0 if y <= LOW else 1 if y >= 0.75 * VDD else None
    return y, got


def main(argv: list[str]) -> int:
    names = argv or list(GATES)
    failures = 0
    for name in names:
        gate = GATES[name]
        header = " ".join(gate.inputs)
        print(f"\n{gate.name}: {gate.logic}  (VDD = {VDD} V)")
        print(f"  {header} | y_sim   logic  ok?")
        print("  " + "-" * (len(header) + 22))
        for inputs in combos(gate):
            y, got = run_case(gate, inputs)
            want = gate.expect(*inputs)
            ok = got == want
            failures += not ok
            print(
                f"  {''.join(map(str, inputs)):>{len(header)}} | {y:6.3f} V"
                f"   {got}      {'OK' if ok else 'FAIL (want %d)' % want}"
            )
    print()
    if failures:
        print(f"FAILED: {failures} row(s) mismatched")
        return 1
    print("all truth tables verified")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
