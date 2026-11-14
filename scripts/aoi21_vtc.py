#!/usr/bin/env python
"""DC sweep of the AOI21 cell: VTC for two input corners.

B=1 (1.8 V), C=0 always: then y = !((A & 1) | 0) = !A, i.e. the AOI21 becomes
an inverter on input A.  Sweeping A from 0 to 1.8 V traces the transfer curve
and shows the switching threshold near VDD/2.

Usage::

    python scripts/aoi21_vtc.py docs/aoi21_vtc.png
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # repo root

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from minispice import NetlistParser, Simulation
from tools.gates import CELLS, MODELS, VDD


def build_netlist(b: float, c: float) -> str:
    return "\n".join(
        [
            "* AOI21 VTC: y = !((a & b) | c), b=1, c=0 -> inverter on a",
            f".include {MODELS}",
            f".include {CELLS / 'aoi21.cir'}",
            f"VDD vdd 0 {VDD}",
            f"VB b 0 {b}",
            f"VC c 0 {c}",
            "VA a 0 0",
            "XG a b c y vdd aoi21",
            ".dc VA 0 1.8 0.02",
            ".end",
        ]
    ) + "\n"


def main(argv: list[str]) -> int:
    out = Path(argv[0]) if argv else Path("docs/aoi21_vtc.png")
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for label, b in (("corner: b = 1 (a-inverter), c = 0", 1.0),):
        circuit = NetlistParser(build_netlist(b, 0.0)).parse()
        sweep = Simulation(circuit).dc_sweep("VA", 0.0, 1.8, 0.02)
        y = sweep.nodes["y"]
        ax.plot(sweep.sweep, y, label=label)
        mid = np.argmin(np.abs(y - VDD / 2))
        ax.plot(sweep.sweep[mid], y[mid], "o", ms=5)
        print(f"switching point: VA = {sweep.sweep[mid]:.2f} V")

    # Ideal inverter reference line
    ax.plot([0, VDD], [VDD, 0], "k--", lw=1, label="ideal !a")
    ax.set_xlabel("VA (V)")
    ax.set_ylabel("V(y)")
    ax.set_title("AOI21 voltage transfer curve (sky130-inspired Level-1)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=140)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
