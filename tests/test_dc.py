import math

import pytest

from minispice import NetlistParser, Simulation


def test_resistive_divider(run):
    result = run(
        "V1 1 0 10\nR1 1 2 1k\nR2 2 0 1k\n.op\n.end\n"
    )[0]
    assert result.nodes["2"] == pytest.approx(5.0, abs=1e-6)
    assert result.currents["V1"] == pytest.approx(-5e-3, rel=1e-6)


def test_current_source(run):
    result = run("I1 0 1 1m\nR1 1 0 1k\n.op\n.end\n")[0]
    assert result.nodes["1"] == pytest.approx(1.0, abs=1e-6)


def test_inductor_is_short_in_dc(run):
    result = run(
        "V1 1 0 1\nR1 1 2 1k\nL1 2 0 1m\n.op\n.end\n"
    )[0]
    assert result.nodes["2"] == pytest.approx(0.0, abs=1e-6)
    assert result.currents["L1"] == pytest.approx(1e-3, rel=1e-6)


def test_dc_sweep(run):
    result = run(
        "V1 1 0 0\nR1 1 2 1k\nR2 2 0 1k\n.dc V1 0 10 2\n.end\n"
    )[0]
    assert list(result.sweep) == pytest.approx([0, 2, 4, 6, 8, 10])
    assert result.nodes["2"] == pytest.approx(result.sweep / 2)


def test_subcircuit_operating_point(run):
    result = run(
        "V1 1 0 10\nX1 1 2 divider\n"
        ".subckt divider in out\nR1 in out 1k\nR2 out 0 1k\n.ends\n"
        ".op\n.end\n"
    )[0]
    assert result.nodes["2"] == pytest.approx(5.0, abs=1e-6)


def test_zero_resistance_is_rejected(run):
    """A 0-ohm resistor has no finite DC stamp (1/R would blow up)."""
    with pytest.raises(ValueError, match="zero resistance"):
        run("V1 1 0 10\nR1 1 2 0\nR2 2 0 1k\n.op\n.end\n")


def test_meas_extrema_over_single_point_dc_sweep(run):
    """MAX/MIN .meas are valid on a degenerate (start == stop) sweep.

    Reference: divider V(out) = V1 * R2/(R1+R2) = 0.5 V exactly for every
    sweep point, so MAX = MIN = 0.5 V with one point.  A blanket
    "at least two points" rule previously turned this into NaN (yield
    benches with .dc VDD 1.8 1.8 1 hit it).  Window kinds (RISE/FALL/
    CROSS/WIDTH) still require >= 2 points; AVG/RMS still reject a
    zero-length interval.
    """
    result = run(
        "V1 in 0 1\nR1 in out 1k\nR2 out 0 1k\n"
        ".dc V1 1 1 1\n"
        ".meas dc vmax MAX V(out) MIN=0.4 MAX=0.6\n"
        ".meas dc vmin MIN V(out)\n"
        ".end\n"
    )[0]
    assert result.measurements["vmax"] == pytest.approx(0.5, abs=1e-9)
    assert result.measurements["vmin"] == pytest.approx(0.5, abs=1e-9)

    # Window kinds still need >= 2 points even though MAX/MIN do not;
    # the meas engine converts the ValueError into NaN + a warning.
    circuit = NetlistParser(
        "V1 in 0 1\nR1 in out 1k\nR2 out 0 1k\n"
        ".dc V1 1 1 1\n"
        ".meas dc trise RISE V(out)\n.end\n"
    ).parse()
    result = Simulation(circuit).run()[0]
    assert math.isnan(result.measurements["trise"])
    assert any("at least two points" in w for w in circuit.warnings)


def test_yield_on_single_point_dc_sweep(run):
    """Monte Carlo yield specs work with a degenerate sweep (regression
    for the zero-populated yield_summary found by the cell verification
    bench on the 34-FET full adder).

    With +/-10% uniform tolerances the divider output spans
    [0.9/2.0, 1.1/1.8] = [0.45, 0.6111] V, so a spec window of
    [0.4, 0.65] V contains every draw: yield 100% (100/100) is the
    deterministic expectation under the fixed seed.
    """
    circuit = NetlistParser(
        "V1 in 0 1\nR1 in out 1k TOL=0.1\nR2 out 0 1k TOL=0.1\n"
        ".dc V1 1 1 1\n"
        ".meas dc vout MAX V(out) MIN=0.4 MAX=0.65\n"
        ".mc 100 SEED=42\n.end\n"
    ).parse()
    sim = Simulation(circuit)
    sim.run()
    summary = sim.yield_summary["vout"]
    assert summary["total"] == 100
    assert summary["pass"] == 100
    assert summary["yield"] == pytest.approx(1.0)
