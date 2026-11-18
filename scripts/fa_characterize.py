#!/usr/bin/env python
"""Characterize sky130 full adder (fa): delays vs CL, leakage, corners.

Produces reproducible numbers for docs/fa_characterization.md.
Run: .venv/bin/python scripts/fa_characterize.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np

from minispice import NetlistParser, Simulation
from tools.gates import CELLS, MODELS, VDD

CL_LIST = [2e-15, 4e-15, 8e-15, 16e-15]
EDGE_PS = 100e-12
TSTEP = 1e-12
TSTOP = 6e-9


def measure_50pct(t, v, vdd, t_search_start=0.3e-9):
    """Return dict with t_r50, t_f50 (first rising/falling 50% after start)."""
    th = 0.5 * vdd
    idx0 = np.searchsorted(t, t_search_start)
    seg_t = t[idx0:]
    seg_v = v[idx0:]
    out = {}
    # rising
    for i in range(1, len(seg_v)):
        if seg_v[i - 1] < th <= seg_v[i]:
            frac = (th - seg_v[i - 1]) / (seg_v[i] - seg_v[i - 1] + 1e-30)
            out["t_r50"] = float(seg_t[i - 1] + frac * (seg_t[i] - seg_t[i - 1]))
            break
    # falling
    for i in range(1, len(seg_v)):
        if seg_v[i - 1] > th >= seg_v[i]:
            frac = (th - seg_v[i - 1]) / (seg_v[i] - seg_v[i - 1] + 1e-30)
            out["t_f50"] = float(seg_t[i - 1] + frac * (seg_t[i] - seg_t[i - 1]))
            break
    return out


def build_delay_netlist(cl_sum, cl_cout, vdd=VDD, temp_c=27.0):
    rise = EDGE_PS
    fall = EDGE_PS
    delay = 0.4e-9
    width = 1.2e-9
    period = 3.0e-9
    lines = [
        "* fa delay: A=0,B=1, CIN pulses 0->1.8, edges 100ps",
        f".include {MODELS}",
        f".include {CELLS / 'fa.cir'}",
        f"VDD vdd 0 {vdd}",
        "VA a 0 0",
        f"VB b 0 {vdd}",
        f"VCIN cin 0 PULSE(0 {vdd} {delay} {rise} {fall} {width} {period})",
        "XG a b cin sum cout vdd fa",
        f"CLSUM sum 0 {cl_sum}",
        f"CLCOUT cout 0 {cl_cout}",
        f".temp {temp_c}",
        f".tran {TSTEP} {TSTOP}",
        ".end",
    ]
    return "\n".join(lines) + "\n"


def run_delays_for_cl(cl, vdd=VDD, temp_c=27.0):
    net = build_delay_netlist(cl, cl, vdd, temp_c)
    c = NetlistParser(net).parse()
    sim = Simulation(c)
    tr = sim.transient(TSTEP, TSTOP)
    t_in50_rise = 0.4e-9 + 0.5 * EDGE_PS
    t_in50_fall = 0.4e-9 + 1.2e-9 + 0.5 * EDGE_PS
    ms = measure_50pct(tr.time, tr.nodes["sum"], vdd, t_search_start=0.45e-9)
    mc = measure_50pct(tr.time, tr.nodes["cout"], vdd, t_search_start=0.45e-9)
    res = {"cl_f": cl}
    if "t_f50" in ms:
        res["sum_tpHL"] = ms["t_f50"] - t_in50_rise
    if "t_r50" in ms:
        res["sum_tpLH"] = ms["t_r50"] - t_in50_fall
    if "t_r50" in mc:
        res["cout_tpLH"] = mc["t_r50"] - t_in50_rise
    if "t_f50" in mc:
        res["cout_tpHL"] = mc["t_f50"] - t_in50_fall
    return res
