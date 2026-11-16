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


def test_capacitor_ic_with_grounded_first_terminal(run):
    """`C1 0 out 1u IC=-0.5` sets v(n1) - v(n2) = IC, so v(out) = +0.5 V at t=0.

    Same RC decay as `test_rc_step` but with the capacitor written the "wrong"
    way round, which is the branch of the UIC initial-state code that negates
    the IC because n1 is ground.
    """
    result = run(
        "C1 0 out 1u IC=-0.5\n"
        "R1 out 0 1k\n"
        ".tran 10u 1m uic\n.end\n"
    )[0]
    assert result.nodes["out"][0] == pytest.approx(0.5, abs=1e-12)
    # tau = RC = 1 ms, so v(1 ms) = 0.5*exp(-1) = 0.183940 V; dt = 10 us gives
    # N = 100 steps and a ~0.3% accumulated truncation error (measured 0.5%).
    assert result.nodes["out"][-1] == pytest.approx(0.5 / np.e, abs=2e-3)


def test_inductor_ic_sets_initial_branch_current(run):
    """`IC=` on an inductor with `.tran ... uic` seeds iL(0) in the MNA branch.

    R1 || L1 with an initial current and no source is a pure RL decay:
    iL(t) = I0*exp(-t/tau), tau = L/R = 10 mH / 10 ohm = 1 ms, I0 = 1 A.
    """
    result = run(
        "L1 out 0 10m IC=1\n"
        "R1 out 0 10\n"
        ".tran 2u 1m uic\n.end\n"
    )[0]
    tau = 10e-3 / 10.0
    # The initial row is the analytic initial state: iL(0) = 1 A, v(out) = 0.
    assert result.time[0] == 0.0
    assert result.currents["L1"][0] == pytest.approx(1.0, abs=1e-12)
    assert result.nodes["out"][0] == pytest.approx(0.0, abs=1e-12)
    # One time constant later iL = exp(-1) = 0.367879 A. With dt = 2 us the
    # local truncation error is O((dt/tau)^2) = 4e-6 per step; over 500 steps
    # that accumulates to ~1e-3 relative, so rel=2e-3 covers the discretization
    # (measured deviation is 1.0e-3).
    at_tau = result.currents["L1"][np.argmin(np.abs(result.time - tau))]
    assert at_tau == pytest.approx(np.exp(-1.0), rel=2e-3)
    final = result.currents["L1"][-1]
    assert final == pytest.approx(np.exp(-1.0), rel=2e-3)
    # Half a time constant: iL = exp(-0.5) = 0.606531 A.
    at_half = result.currents["L1"][np.argmin(np.abs(result.time - 0.5 * tau))]
    assert at_half == pytest.approx(np.exp(-0.5), rel=2e-3)
    # KCL at `out`: the same current flows through R1, and L1 is wired out -> 0,
    # so v(out) = -R*I(L1) = -3.68248 V exactly (identity, not an approximation).
    assert result.nodes["out"][-1] == pytest.approx(-10.0 * final, abs=1e-9)


def test_transient_halves_step_after_convergence_failure(run):
    """A failing Newton step is retried at dt/2 and must still reach the answer.

    `.options maxiter=2` makes the Newton solve genuinely fail on the stiff
    diode step (the capacitor node jumps by hundreds of mV per 100 us step), so
    the transient loop catches ConvergenceError and halves dt before retrying.
    The final value must still be the analytic steady state: 1 mA forced through
    D1 by I1 (no current left for C1), so V = Vt*ln(1 + I/IS).
    """
    result = run(
        "Diode clamp\n"
        "I1 0 out PULSE(0 1m 0 1n 1n 10m 20m)\n"
        "C1 out 0 1u\n"
        "D1 out 0 dmod\n"
        ".model dmod D(IS=1e-12 N=1)\n"
        ".options maxiter=2\n"
        ".tran 100u 4m uic\n.end\n"
    )[0]
    steps = np.diff(result.time)
    # Evidence that the retry branch ran: the loop restarts at dt/2 (<= 50 us),
    # so the grid is finer than half the 100 us maximum step and contains more
    # rows than the 4 ms / 100 us = 40 uniform steps would give.
    assert steps.min() < 0.5 * steps.max()
    assert len(result.time) > 41
    # Vt = kT/q at the default 27 C = 1.380649e-23 * 300.15 / 1.602176634e-19
    # = 25.8649 mV; V = Vt*ln(1 + 1 mA / 1e-12 A) = 0.536006 V.
    vt = 1.380649e-23 * (27.0 + 273.15) / 1.602176634e-19
    analytic = vt * np.log(1.0 + 1e-3 / 1e-12)
    assert result.nodes["out"][-1] == pytest.approx(analytic, rel=1e-5)


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
