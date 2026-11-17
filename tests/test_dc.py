import pytest


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
