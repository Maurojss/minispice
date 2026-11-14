#!/usr/bin/env python3
"""Parasitic capacitance of the sky130 transistors used by our cells.

Follows the Weste & Harris (CMOS VLSI Design, ch. 4 - Circuit
Characterization) methodology, with parameters taken from the real sky130
BSIM4 PM3 files vendored in `models/`:

  1. Gate capacitance:        CGG = Cox*W*L + (cgso+cgdo)*W
  2. Junction capacitance:    Cj(V) = CJ*AD*(1+V/PB)^-MJ
                                       + CJSW*PD*(1+V/PBSW)^-MJSW
                                       + CJSWG*PG*(1+V/PBSWG)^-MJSWG
                              (grade-coefficient form; evaluated at VDD/2)
  3. Effective resistance:    Req = (1/0.69) * int(VDD/2 .. VDD) dv / Id(v)
                              (saturation-current / alpha-power framework)
  4. Delay estimate:          t = 0.69 * Req * CL   -> verified in-sim.

Usage:
    python pdk/sky130/capacitance.py [--D 0.5] [--VDD 1.8]
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

EPS_0 = 8.854187817e-12  # F/m
EPS_OX = 3.9 * EPS_0
MODELS = Path(__file__).resolve().parent / "models"

# Diffusion height assumed for the S/D strip (documented teaching value:
# a contacted diffusion tap in a cell layout; rescale for your layout).
D_DIFF = 0.5e-6
VDD = 1.8


def grep_first(text: str, param: str) -> float:
    m = re.search(rf"\+ {param} = ([\d.e+\-]+)", text)
    if not m:
        raise ValueError(f"parameter {param} not found")
    return float(m.group(1))


@dataclass
class DeviceCaps:
    """Capacitance-relevant parameters of one sky130 device (tt corner)."""

    label: str
    w: float
    l: float
    toxm: float
    cgso: float  # F/m (per width), source & drain side overlap
    cgdo: float
    cjs: float  # F/m^2, zero-bias area junction
    mjs: float
    pbs: float
    cjsws: float  # F/m, lateral sidewall
    mjsws: float
    pbsws: float
    cjswgs: float  # F/m, gate-edge sidewall
    mjswgs: float
    pbswgs: float

    # ----------------------------------------------------------- Level-1 R
    kp: float
    vto: float
    lam: float

    # ------------------------------------------------------------- caps
    @property
    def cox(self) -> float:
        return EPS_OX / self.toxm

    def gate_cap(self) -> float:
        """Total gate capacitance in F (channel + both overlaps)."""
        return self.cox * self.w * self.l + (self.cgso + self.cgdo) * self.w

    def gate_cap_per_micron(self) -> float:
        return self.gate_cap() / self.w * 1e-15 / 1e-6

    def junction_cap(self, v: float, d: float = D_DIFF) -> float:
        """Drain-junction capacitance at reverse bias v (Harris grading form).

        Geometry: diffusion strip of width W, height d.
          AD  = W*d          (area term)
          PD  = 2d + W       (three outer edges: two sides + far edge)
          PG  = W            (edge facing the gate: gate-edge sidewall)
        """
        ad = self.w * d
        pd = 2.0 * d + self.w
        pg = self.w
        return (
            self.cjs * ad / (1.0 + v / self.pbs) ** self.mjs
            + self.cjsws * pd / (1.0 + v / self.pbsws) ** self.mjsws
            + self.cjswgs * pg / (1.0 + v / self.pbswgs) ** self.mjswgs
        )

    # ------------------------------------------------------ 0.69*Req*C
    def req_eff(self, vdd: float) -> float:
        """Effective pull-down/pull-up resistance of Level-1 model.

        Req = (1/0.69) * integral_{vdd/2}^{vdd} dv / Id(vds=v),
        Id(v) from the simulator's own MOS model (sat + linear + lambda).
        """
        beta = self.kp * (self.w / self.l)
        vth = abs(self.vto)
        sign = 1.0 if self.vto > 0 else -1.0  # pfet: pick pfet sign internally

        def id_at(vds: float) -> float:
            sign_i = 1.0 if self.kp > 0 else -1.0
            vov = vdd - vth
            if vov <= 0:
                return 1e-9
            if vds >= vov:
                return 0.5 * abs(beta) * vov * vov * (1.0 + self.lam * vds)
            return abs(beta) * (vov * vds - 0.5 * vds * vds) * (
                1.0 + self.lam * vds
            )

        n = 2000
        total = 0.0
        for i in range(n):
            a = vdd / 2.0 + (vdd / 2.0) * i / n
            b = vdd / 2.0 + (vdd / 2.0) * (i + 1) / n
            total += 0.5 * (b - a) * (1.0 / id_at(a) + 1.0 / id_at(b))
        return total / 0.69


def load(label: str, fname: str, w: float, l: float, kp: float, vto: float, lam: float) -> DeviceCaps:
    text = (MODELS / fname).read_text(encoding="utf-8", errors="ignore")
    return DeviceCaps(
        label=label,
        w=w,
        l=l,
        toxm=grep_first(text, "toxm"),
        cgso=grep_first(text, "cgso"),
        cgdo=grep_first(text, "cgdo"),
        cjs=grep_first(text, "cjs"),
        mjs=grep_first(text, "mjs"),
        pbs=grep_first(text, "pbs"),
        cjsws=grep_first(text, "cjsws"),
        mjsws=grep_first(text, "mjsws"),
        pbsws=grep_first(text, "pbsws"),
        cjswgs=grep_first(text, "cjswgs"),
        mjswgs=grep_first(text, "mjswgs"),
        pbswgs=grep_first(text, "pbswgs"),
        kp=kp,
        vto=vto,
        lam=lam,
    )


def fmt(f: float) -> str:
    return f"{f * 1e15:6.3f}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--D", type=float, default=D_DIFF, help="diffusion height (m)")
    ap.add_argument("--vdd", type=float, default=VDD)
    args = ap.parse_args()
    d, vdd = args.D, args.vdd

    nfet = load("NFET (W=1.0 um, L=0.15 um)", "nfet_tt.pm3.spice", 1.0e-6, 0.15e-6,
                kp=240e-6, vto=0.52, lam=0.03)
    pfet = load("PFET (W=3.0 um, L=0.15 um)", "pfet_tt.pm3.spice", 3.0e-6, 0.15e-6,
                kp=-86e-6, vto=-0.75, lam=0.05)

    print("sky130 parasitic capacitance (Weste & Harris method, tt corner)")
    print(f"tox: nfet {nfet.toxm*1e9:.3f} nm, pfet {pfet.toxm*1e9:.3f} nm")
    print(f"Cox: nfet {nfet.cox*1e3:.3f} mF/m^2, pfet {pfet.cox*1e3:.3f} mF/m^2")
    print(f"diffusion height D = {d*1e6:.2f} um, VDD = {vdd} V, caps at V = VDD/2")
    print()
    print(f"{'device':26s} {'CGG [fF]':>9s} {'Cj(VDD/2)':>10s} {'Cj(0)':>8s}")
    for dev in (nfet, pfet):
        print(
            f"{dev.label:26s} {fmt(dev.gate_cap()):>9s} "
            f"{fmt(dev.junction_cap(vdd / 2, d)):>10s} "
            f"{fmt(dev.junction_cap(0.0, d)):>8s}"
        )
    print()

    cg_nf, cg_pf = nfet.gate_cap(), pfet.gate_cap()
    cj_nf = nfet.junction_cap(vdd / 2, d)
    cj_pf = pfet.junction_cap(vdd / 2, d)

    # ------- load of one fanout inverter (the compound cells are FO=1 chains)
    fov = cg_nf + cg_pf  # gate cap of the following inverter
    cout = cj_nf + cj_pf  # self-loading at the output node
    cl = cout + fov
    print("inverter (FO1, no wire load):")
    print(f"  gate load of next inverter : {fmt(fov)} fF")
    print(f"  self-loading (S/D junction): {fmt(cout)} fF")
    print(f"  total CL                   : {fmt(cl)} fF")
    print()

    rn = nfet.req_eff(vdd)
    rp = pfet.req_eff(vdd)
    tphl = 0.69 * rn * cl
    tplh = 0.69 * rp * cl
    print(f"Req (0.69*Req*CL basis):  NMOS {rn:7.0f} ohm, PMOS {rp:7.0f} ohm")
    print(f"predicted delays:  tPHL = {tphl*1e12:5.2f} ps, tPLH = {tplh*1e12:5.2f} ps")

    # -------- netlist-ready equivalent capacitances for the transient example
    print()
    print("netlist-ready values (examples/cmos_inverter_cap.cir):")
    print(f"  Cg_n = {cg_nf*1e15:.2f}f  Cg_p = {cg_pf*1e15:.2f}f  (next-stage gate caps)")
    print(f"  Cout ≈ {(cout+fov)*1e15:.2f}f at the output node")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
