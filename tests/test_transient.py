import numpy as np
import pytest


def test_rc_step(run):
    result = run(
        "V1 in 0 PULSE(0 1 0 1n 1n 100m 200m)\n"
        "R1 in out 1k\nC1 out 0 1u\n"
        ".tran 10u 5m\n.end\n"
    )[0]
    final = result.nodes["out"][-1]
    at_tau = result.nodes["out"][np.argmin(np.abs(result.time - 1e-3))]
    assert final == pytest.approx(1.0 - np.exp(-5), abs=5e-3)
    assert at_tau == pytest.approx(1.0 - np.exp(-1), abs=2e-2)


def test_rl_step(run):
    result = run(
        "V1 in 0 PULSE(0 1 0 1n 1n 100m 200m)\n"
        "R1 in out 10\nL1 out 0 1m\n"
        ".tran 2u 1m\n.end\n"
    )[0]
    final = result.currents["L1"][-1]
    at_tau = result.currents["L1"][np.argmin(np.abs(result.time - 0.1e-3))]
    assert final == pytest.approx(0.1, abs=2e-3)
    assert at_tau == pytest.approx(0.1 * (1 - np.exp(-1)), abs=3e-3)


def test_rlc_overshoot(run):
    result = run(
        "V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n"
        "R1 in n1 10\nL1 n1 out 1m\nC1 out 0 1u\n"
        ".tran 1u 2m\n.end\n"
    )[0]
    assert np.max(result.nodes["out"]) > 1.2
    assert result.nodes["out"][-1] == pytest.approx(1.0, abs=1e-2)


def test_backward_euler_option(run):
    result = run(
        "V1 in 0 PULSE(0 1 0 1n 1n 100m 200m)\n"
        "R1 in out 1k\nC1 out 0 1u\n"
        ".options method=be\n"
        ".tran 10u 5m\n.end\n"
    )[0]
    assert result.nodes["out"][-1] == pytest.approx(1.0 - (1 / 1.1) ** 50, abs=2e-3)


def test_cmos_inverter_transient(run):
    result = run(
        "VDD vdd 0 1.8\n"
        "VIN in 0 PULSE(0 1.8 0 100p 100p 2n 4n)\n"
        "MP out in vdd vdd PMOS W=2u L=1u\n"
        "MN out in 0 0 NMOS W=1u L=1u\n"
        "CL out 0 10f\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
        ".model PMOS PMOS(VTO=-0.5 KP=60u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
        ".tran 10p 8n\n.end\n"
    )[0]
    at_1ns = result.nodes["out"][np.argmin(np.abs(result.time - 1e-9))]
    at_3ns = result.nodes["out"][np.argmin(np.abs(result.time - 3e-9))]
    assert at_1ns < 0.1
    assert at_3ns > 1.7
