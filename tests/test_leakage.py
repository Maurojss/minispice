"""Leakage-current tests: diode reverse saturation + MOSFET subthreshold."""

import math

import numpy as np

from minispice import NetlistParser
from minispice.solver import Simulation

# gmin so pequeno para nao contaminar as medidas de fuga isoladas
_OPTIONS = ".options gmin=1e-15"


def test_diode_reverse_leakage_matches_isat():
    """Off-state diode current exposed per element: I(D1) ~ -IS."""
    netlist = f"""diode leakage
V1 a 0 dc -1
D1 a b DMOD
R1 b 0 1meg
.model DMOD D(IS=1e-12 N=1.0)
{_OPTIONS}
.op
.end
"""
    op = Simulation(NetlistParser(text=netlist).parse()).operating_point()
    assert "D1" in op.currents                     # exposed on its own
    assert abs(op.currents["D1"] - (-1e-12)) < 1e-15   # saturates at -IS
    assert abs(op.currents["V1"] - 1e-12) < 5e-15  # source carries the leak


def test_diode_forward_current_reported():
    netlist = f"""diode forward
V1 a 0 dc 1
D1 a b DMOD
R1 b 0 1k
.model DMOD D(IS=1e-12)
{_OPTIONS}
.op
.end
"""
    op = Simulation(NetlistParser(text=netlist).parse()).operating_point()
    # IS=1e-12, n=1 gives a ~0.52 V drop: I = (1 - v_diode)/1k ~ 0.48 mA
    i = op.currents["D1"]
    assert 0.3e-3 < abs(i) < 0.7e-3
    assert abs(abs(op.nodes["b"]) - abs(i) * 1e3) < 1e-7   # V(b) = I * R1
    # KCL stays consistent to Newton's current tolerance (reltol * |I|)
    assert abs(i + op.currents["V1"]) < abs(i) * 1e-6


def test_diode_leakage_sweep_is_visible_as_i_label():
    """The GUI/CLI series plumbing: dc_sweep carries I:D1 automatically."""
    netlist = f"""diode leakage sweep
V1 a 0 dc -0.5
D1 a b D1MOD
R1 b 0 1meg
.model D1MOD D(IS=2e-12)
{_OPTIONS}
.dc V1 -2 -0.5 0.5
.end
"""
    sweep = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    assert "D1" in sweep.currents
    series = np.asarray(sweep.currents["D1"])
    assert len(series) == 4
    # every point in reverse bias saturates at +IS flowing out of the anode
    assert np.allclose(series, 2e-12, rtol=1e-6)


def test_mosfet_subthreshold_leakage():
    """Off-state FET leaks: Id = 0.5*beta*(nsub*vt)^2*exp(vov/(nsub*vt))."""
    vt = 1.380649e-23 * 300.15 / 1.602176634e-19       # 27 C
    vp = 1.5 * vt
    beta = 120e-6
    netlist = f"""mosfet leakage
VDD drain 0 1.8
Vg gate 0 0
MN drain gate 0 0 NMOS W=1u L=1u
.model NMOS NMOS(VTO=0.5 KP=120u)
{_OPTIONS}
.op
.end
"""
    circuit = NetlistParser(text=netlist).parse()
    op = Simulation(circuit).operating_point()
    expected = 0.5 * beta * vp**2 * math.exp(-0.5 / vp)
    i = abs(op.currents["VDD"])
    # gmin (~1.8e-15) rides along on the ammeter branch
    assert abs(i - expected) < 2e-15
    assert 1e-14 < i < 1e-12                          # femto/pico scale


def test_subthreshold_slope_is_one_decade_per_nsub_vt_gene():
    """Id(x) / Id(x - nsub*vt*ln(10)) == 10 exactly on the model's own grid."""
    vt = 1.380649e-23 * 300.15 / 1.602176634e-19
    base = """mosfet subth
VDD drain 0 1.8
Vg gate 0 {vgs}
MN drain gate 0 0 NMOS W=1u L=1u
.model NMOS NMOS(VTO=0.5 KP=120u NSUB=1.5)
.op
.end
"""
    def leak(vgs: float) -> float:
        netlist = base.replace("{vgs}", f"{vgs:.12g}") + _OPTIONS
        op = Simulation(NetlistParser(text=netlist).parse()).operating_point()
        return abs(op.currents["VDD"])

    step = 1.5 * vt * math.log(10)
    ratio = leak(0.5) / leak(0.5 - step)
    assert abs(ratio - 10.0) < 1e-3


def test_nsub_zero_disables_leakage():
    """NSUB=0 keeps the classic Level-1 hard cutoff (no leak)."""
    base = """mosfet nsub
VDD drain 0 1.8
Vg gate 0 0
MN drain gate 0 0 NMOS W=1u L=1u
.model NMOS NMOS(VTO=0.5 KP=120u {nsub})
.op
.end
"""
    def leak(nsub_line: str) -> float:
        netlist = base.replace("{nsub}", nsub_line) + _OPTIONS
        op = Simulation(NetlistParser(text=netlist).parse()).operating_point()
        return abs(op.currents["VDD"])

    with_leak = leak("NSUB=1.5")
    without = leak("NSUB=0")
    assert with_leak > 1e-13
    assert without < 1e-14                          # only the tiny gmin path


def test_cmos_pmos_leakage_scales_with_width():
    """Off-state PMOS leak: I0 scales with W (3x width -> 3x leak)."""
    vt = 1.380649e-23 * 300.15 / 1.602176634e-19
    vp = 1.5 * vt
    i0 = 0.5 * 60e-6 * vp**2

    def leak(width: str) -> float:
        netlist = f"""pmos leak
VDD vdd 0 1.8
VIN in 0 1.8
MP out in vdd vdd PMOS {width} L=1u
R1 out 0 1meg
.model PMOS PMOS(VTO=-0.5 KP=60u)
{_OPTIONS}
.op
.end
"""
        op = Simulation(NetlistParser(text=netlist).parse()).operating_point()
        return abs(op.currents["VDD"])

    narrow = leak("W=1u")     # I0 * exp(-0.5/vp)
    wide = leak("W=3u")       # 3x I0 * exp(-0.5/vp)
    assert narrow > i0 * math.exp(-0.5 / vp) * 0.5
    assert abs(wide - 3 * narrow) < 3 * 5e-15          # width ratio preserved


def test_gui_signals_include_diode_currents():
    """What the GUI receives: I:D1 sits next to I:V1 in the signal list."""
    from minispice.webui import _result_payload

    netlist = f"""gui check
V1 a 0 dc -1
D1 a b DMOD
R1 b 0 1meg
.model DMOD D(IS=1e-12)
{_OPTIONS}
.dc V1 -2 -2 1
.end
"""
    result = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    payload = _result_payload(0, result)
    labels = [sig["label"] for sig in payload["signals"]]
    assert "I:D1" in labels and "I:V1" in labels and "V:a" in labels
