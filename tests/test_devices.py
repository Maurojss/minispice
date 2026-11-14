import numpy as np
import pytest


def test_diode_operating_point(run):
    result = run(
        "V1 in 0 1\nR1 in out 1k\nD1 out 0 DMOD\n"
        ".model DMOD D(IS=1e-14 N=1)\n.op\n.end\n"
    )[0]
    assert 0.55 < result.nodes["out"] < 0.72
    current = (1.0 - result.nodes["out"]) / 1e3
    assert result.currents["V1"] == pytest.approx(-current, rel=1e-6)


def test_diode_rectifier(run):
    result = run(
        "V1 in 0 SIN(0 5 1k)\nD1 in out DMOD\nR1 out 0 1k\n"
        ".model DMOD D(IS=1e-14 N=1)\n.tran 10u 3m\n.end\n"
    )[0]
    assert np.max(result.nodes["out"]) == pytest.approx(4.3, abs=0.3)
    assert np.min(result.nodes["out"]) > -0.05


def test_cmos_inverter_vtc(run):
    result = run(
        "VDD vdd 0 1.8\nVIN in 0 0\n"
        "MP out in vdd vdd PMOS W=2u L=1u\n"
        "MN out in 0 0 NMOS W=1u L=1u\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
        ".model PMOS PMOS(VTO=-0.5 KP=60u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
        ".dc VIN 0 1.8 0.01\n.end\n"
    )[0]
    out = result.nodes["out"]
    assert out[0] > 1.75
    assert out[-1] < 0.05
    middle = out[np.argmin(np.abs(result.sweep - 0.9))]
    assert 0.7 < middle < 1.1


def test_diode_reverse_blocking(run):
    result = run(
        "V1 in 0 -5\nR1 in out 1k\nD1 out 0 DMOD\n"
        ".model DMOD D(IS=1e-14 N=1)\n.op\n.end\n"
    )[0]
    assert result.nodes["out"] == pytest.approx(-5.0, abs=0.01)
    assert abs(result.currents["V1"]) < 1e-9
