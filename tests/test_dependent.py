"""MNA stamping of dependent and behavioral sources (E/G/F/H/B) in DC and AC.

Every circuit here has a closed-form answer derived by hand; the numbers in the
assertions cite that derivation. The DC half exercises the differential
stamping forms (both control terminals and both output terminals away from
ground) plus the explicit current branch that E/H carry, and the AC half checks
that the same stamps give the hand-derived complex transfer.
"""

import math

import pytest

from minispice import NetlistParser
from minispice.solver import Simulation

# Frequency at which omega*R*C = 1000 rad/s * 1 kohm * 1 uF = 1 exactly, so the
# RC impedances below become the pure imaginary unit j and the algebra stays
# hand-checkable.
_OMEGA_RC = 1000.0
_F_RC = _OMEGA_RC / (2.0 * math.pi)


def _run(text: str):
    return Simulation(NetlistParser(text).parse()).run()


def test_vcvs_differential_control_and_branch_current_dc():
    """E1 with both control pins off ground, and its current branch observable.

    E1 forces V(out) - V(ref) = 2*(V(in) - V(ref)). Neither control terminal is
    ground, so the stamp must write both A[k, nc1] and A[k, nc2].
    """
    result = _run(
        "VCVS differential\n"
        "V1 in 0 5\n"
        "V2 ref 0 1\n"
        "E1 out ref in ref 2\n"
        "R1 out 0 1k\n"
        "R2 ref 0 1k\n"
        ".op\n.end\n"
    )[0]
    # Closed form: V(out) = V(ref) + 2*(V(in) - V(ref)) = 1 + 2*(5 - 1) = 9 V.
    assert result.nodes["out"] == pytest.approx(9.0, abs=1e-6)
    assert result.nodes["ref"] == pytest.approx(1.0, abs=1e-6)
    # KCL at `out`: only R1 and the E1 branch meet there, so the branch current
    # (defined n1 -> n2, i.e. `out` -> `ref` inside the source) cancels R1:
    #   I_E1 = -V(out)/1k = -9 mA.
    assert result.currents["E1"] == pytest.approx(-9e-3, abs=1e-9)
    # Consistency check at `ref`: R2 sinks 1 mA, E1 pushes 9 mA back in, so the
    # stiff V2 branch must carry I_V2 = -(9 mA + 1 mA) = -10 mA.
    assert result.currents["V2"] == pytest.approx(-10e-3, abs=1e-9)


def test_vccs_floating_output_and_differential_control_dc():
    """G1 delivers gm*(V(a) - V(b)) between two nodes that are both off ground.

    This is the stamp that touches all four rows/columns at once: both output
    nodes and both control nodes are non-ground.
    """
    result = _run(
        "VCCS floating\n"
        "V1 a 0 2\n"
        "V2 b 0 1\n"
        "G1 x y a b 1m\n"
        "R1 x 0 1k\n"
        "R2 y 0 1k\n"
        ".op\n.end\n"
    )[0]
    # Control voltage V(a) - V(b) = 1 V -> G1 sources 1 mA from x to y.
    # KCL at x: V(x)/1k + 1 mA = 0 -> V(x) = -1 V.
    # KCL at y: V(y)/1k - 1 mA = 0 -> V(y) = +1 V.
    # gmin = 1e-12 S shunts each node through 1 kohm, i.e. a 1e-9 V offset;
    # abs=1e-6 is 1000x that, so it only guards against stamping errors.
    assert result.nodes["x"] == pytest.approx(-1.0, abs=1e-6)
    assert result.nodes["y"] == pytest.approx(1.0, abs=1e-6)


def test_cccs_floating_return_scales_control_branch_current_dc():
    """F1 mirrors k*I(V1) into a load that is not grounded on either side."""
    result = _run(
        "CCCS floating\n"
        "V1 a 0 1\n"
        "R1 a 0 1k\n"
        "F1 x y V1 2\n"
        "R2 x 0 1k\n"
        "R3 y 0 1k\n"
        ".op\n.end\n"
    )[0]
    # V1 drives 1 mA into R1, so the current inside V1 (a -> 0) is I(V1) = -1 mA.
    assert result.currents["V1"] == pytest.approx(-1e-3, abs=1e-9)
    # F1 pushes k*I(V1) = -2 mA from x to y.
    # KCL at x: V(x)/1k - 2 mA = 0 -> V(x) = +2 V.
    # KCL at y: V(y)/1k + 2 mA = 0 -> V(y) = -2 V.
    assert result.nodes["x"] == pytest.approx(2.0, abs=1e-6)
    assert result.nodes["y"] == pytest.approx(-2.0, abs=1e-6)


def test_ccvs_floating_output_scales_control_branch_current_dc():
    """H1 enforces V(x) - V(y) = r*I(V1) with both output pins off ground."""
    result = _run(
        "CCVS floating\n"
        "V1 a 0 1\n"
        "R1 a 0 1k\n"
        "H1 x y V1 1k\n"
        "R2 x 0 1k\n"
        "R3 y 0 1k\n"
        ".op\n.end\n"
    )[0]
    assert result.currents["V1"] == pytest.approx(-1e-3, abs=1e-9)
    # V(x) - V(y) = 1 kohm * I(V1) = -1 V, and the two equal 1 kohm loads split
    # it symmetrically about ground: V(x) = -0.5 V, V(y) = +0.5 V.
    assert result.nodes["x"] == pytest.approx(-0.5, abs=1e-6)
    assert result.nodes["y"] == pytest.approx(0.5, abs=1e-6)
    assert result.nodes["x"] - result.nodes["y"] == pytest.approx(-1.0, abs=1e-6)


def test_vcvs_ac_transfer_across_capacitor_matches_hand_analysis():
    """VCVS in the AC path: E1 senses the voltage across C1 and floats on R3/R4.

    With R = 1 kohm, C = 1 uF and omega = 1000 rad/s, u = j*omega*R*C = j.
    Hand analysis (node a = control +, node b = control -):
      (Va - 1) + u(Va - Vb) = 0            [R1 into `in`]
      u(Vb - Va) + Vb = 0                  [C1 and R2]
      -> Vb = u*Va/(1+u) and Va = (1+u)/(1+2u)
    E1 gives V(out) - V(ref) = 2(Va - Vb); summing the `out` and `ref` KCL rows
    kills the branch current and leaves V(ref)/R4 = 0, so V(ref) = 0 and
    V(out) = 2(Va - Vb).
    """
    result = _run(
        "AC VCVS\n"
        "V1 in 0 dc 0 ac 1 0\n"
        "R1 in a 1k\n"
        "C1 a b 1u\n"
        "R2 b 0 1k\n"
        "E1 out ref a b 2\n"
        "R3 out ref 1k\n"
        "R4 ref 0 1k\n"
        f".ac lin 1 {_F_RC:.12f} {_F_RC:.12f}\n"
        ".end\n"
    )[0]
    assert result.freqs[0] == pytest.approx(_F_RC, rel=1e-9)
    # u = j: Va = (1+j)/(1+2j) = (3-j)/5, Vb = u*Va/(1+u) = (2+j)/5.
    assert result.nodes["a"][0] == pytest.approx(complex(0.6, -0.2), rel=1e-6, abs=1e-9)
    assert result.nodes["b"][0] == pytest.approx(complex(0.4, 0.2), rel=1e-6, abs=1e-9)
    # Va - Vb = (1-2j)/5 -> V(out) = 2(Va - Vb) = 0.4 - 0.8j.
    assert result.nodes["out"][0] == pytest.approx(complex(0.4, -0.8), rel=1e-6, abs=1e-9)
    # V(ref) = 0 by the KCL argument above; the 1e-9 slack is the gmin shunt.
    assert result.nodes["ref"][0] == pytest.approx(0.0 + 0.0j, abs=1e-6)
    # I_E1 = -V(out)/R3 = (-0.4 + 0.8j) mA.
    assert result.currents["E1"][0] == pytest.approx(
        complex(-0.4e-3, 0.8e-3), rel=1e-6, abs=1e-12
    )


def test_ccvs_and_cccs_ac_follow_complex_control_branch_current():
    """H1/F1 in the AC path, controlled by the (complex) current through V1.

    R1 || C1 loads V1, so I(V1) = -(1/R1 + j*omega*C1) = -(1 + j) mA at
    omega*R*C = 1. H1 forces V(hx) - V(hy) = 1 kohm * I(V1) = -(1 + j) V, and
    the equal 1 kohm loads split it as +-(0.5 + 0.5j). F1 injects k*I(V1) =
    -(2 + 2j) mA from fx to fy, so V(fx) = +(2 + 2j) V and V(fy) = -(2 + 2j) V.
    """
    result = _run(
        "AC CCVS CCCS\n"
        "V1 in 0 dc 0 ac 1 0\n"
        "R1 in 0 1k\n"
        "C1 in 0 1u\n"
        "H1 hx hy V1 1k\n"
        "R2 hx 0 1k\n"
        "R3 hy 0 1k\n"
        "F1 fx fy V1 2\n"
        "R4 fx 0 1k\n"
        "R5 fy 0 1k\n"
        f".ac lin 1 {_F_RC:.12f} {_F_RC:.12f}\n"
        ".end\n"
    )[0]
    # I(V1) = -(1 + j) mA: 1 mA through R1 plus 1 mA through C1 at omega*C = 1 mS.
    assert result.currents["V1"][0] == pytest.approx(
        complex(-1e-3, -1e-3), rel=1e-6, abs=1e-12
    )
    assert result.nodes["hx"][0] == pytest.approx(complex(-0.5, -0.5), rel=1e-6, abs=1e-9)
    assert result.nodes["hy"][0] == pytest.approx(complex(0.5, 0.5), rel=1e-6, abs=1e-9)
    assert result.nodes["fx"][0] == pytest.approx(complex(2.0, 2.0), rel=1e-6, abs=1e-9)
    assert result.nodes["fy"][0] == pytest.approx(complex(-2.0, -2.0), rel=1e-6, abs=1e-9)
    # KCL at hx: the H1 branch current is -V(hx)/1k = (0.5 + 0.5j) mA.
    assert result.currents["H1"][0] == pytest.approx(
        complex(0.5e-3, 0.5e-3), rel=1e-6, abs=1e-12
    )


def test_vccs_ac_transconductance_into_rc_load():
    """G1 in the AC path with a complex load: V(x) = gm*Vab/(1/R + j*omega*C).

    The control is a resistive divider, Vab = V(a) - V(b) = 1 - 0.5 = 0.5 V, so
    the source current is 0.5 mA. At omega*R*C = 1 the load admittance is
    1 mS + j*1 mS = 1 mS*(1 + j), hence V(x) = 0.5/(1+j) = 0.25*(1 - j) V, and
    the returning node y carries the same 0.5 mA into 1 kohm -> V(y) = -0.5 V.
    """
    result = _run(
        "AC VCCS\n"
        "V1 a 0 dc 0 ac 1 0\n"
        "R2 a b 1k\n"
        "R3 b 0 1k\n"
        "G1 y x a b 1m\n"
        "R1 x 0 1k\n"
        "Ry y 0 1k\n"
        "C1 x 0 1u\n"
        f".ac lin 1 {_F_RC:.12f} {_F_RC:.12f}\n"
        ".end\n"
    )[0]
    assert result.nodes["b"][0] == pytest.approx(0.5 + 0.0j, rel=1e-6, abs=1e-9)
    assert result.nodes["x"][0] == pytest.approx(complex(0.25, -0.25), rel=1e-6, abs=1e-9)
    assert result.nodes["y"][0] == pytest.approx(-0.5 + 0.0j, rel=1e-6, abs=1e-9)
    # V1 only feeds the 2 kohm divider: I(V1) = -0.5 mA.
    assert result.currents["V1"][0] == pytest.approx(-0.5e-3 + 0.0j, rel=1e-6, abs=1e-12)


def test_behavioral_voltage_source_ac_across_two_non_ground_nodes():
    """B1 sits between `out` and `ref` (neither is ground) with V={2*V(in)}.

    KCL at out: V(out)/1k + I_B = 0; at ref: V(ref)/1k - I_B = 0. Adding them
    gives V(out) + V(ref) = 0, and the source constraint V(out) - V(ref) = 2 V
    (1 V phasor times gain 2) splits it as V(out) = +1 V, V(ref) = -1 V, with
    I_B1 = -V(out)/1k = -1 mA flowing out -> ref.
    Tolerance: the AC solve is exact up to the gmin shunt (1e-12 S against
    1 mS = 1e-9 relative) and roundoff; rel=1e-6 / abs=1e-9 (measured 0.0 V on
    the node phasors and 1.0e-12 A on the branch current).
    """
    result = _run(
        "Behavioral AC\n"
        "V1 in 0 dc 0 ac 1 0\n"
        "B1 out ref V={2*V(in)}\n"
        "R1 out 0 1k\n"
        "R2 ref 0 1k\n"
        ".ac lin 1 1000 1000\n"
        ".end\n"
    )[0]
    assert result.nodes["out"][0] == pytest.approx(1.0 + 0.0j, rel=1e-6, abs=1e-9)
    assert result.nodes["ref"][0] == pytest.approx(-1.0 + 0.0j, rel=1e-6, abs=1e-9)
    assert result.currents["B1"][0] == pytest.approx(-1e-3 + 0.0j, rel=1e-6, abs=1e-12)


@pytest.mark.parametrize("bias", [1.0, 2.5])
def test_behavioral_nonlinear_voltage_source_ac_gain_is_twice_the_bias(bias):
    """V={V(in)**2} linearized around V(in) = VQ gives a small-signal gain 2*VQ.

    d(V(in)^2)/dV(in) = 2*V(in) = 2*VQ at the DC bias, so a 0.1 V AC phasor on
    `in` must appear as 0.2*VQ V on `out` (the B source is stiff, R1 only
    supplies the load current). The .op value itself is VQ^2. The bias is
    parametrized (1 V and 2.5 V) so that a gradient taken at zero instead of at
    the operating point would fail the second case (gain 5 instead of 2).

    Tolerance: the gradient is a central difference with step
    1e-6*max(1,|V|); for a quadratic the truncation error is exactly zero
    (third derivative is 0), leaving float roundoff, so rel=1e-6 / abs=1e-9
    (measured 2.0e-13 V on the phasor and 0.0 V on the DC point).
    """
    results = _run(
        "Behavioral nonlinear\n"
        f"V1 in 0 dc {bias} ac 0.1 0\n"
        "B1 out 0 V={V(in)**2}\n"
        "R1 out 0 1k\n"
        ".op\n"
        ".ac lin 1 1000 1000\n"
        ".end\n"
    )
    op, ac = results[0], results[1]
    assert op.nodes["out"] == pytest.approx(bias ** 2, rel=1e-6, abs=1e-9)
    assert ac.nodes["out"][0] == pytest.approx(2.0 * bias * 0.1 + 0.0j, rel=1e-6, abs=1e-9)


def test_behavioral_current_source_floating_output_dc_and_ac():
    """B2 pushes V(in)**2/1k from p to n, both off ground (no return to 0).

    DC: V(in) = 1 V -> 1 mA from p to n. KCL at p: V(p)/1k + 1 mA = 0 ->
    V(p) = -1 V; at n: V(n)/1k - 1 mA = 0 -> V(n) = +1 V. Both rows need the
    b[j]/A[j,:] stamps, which is what this circuit adds.
    AC: the small-signal transconductance is d/dV(in) = 2*V(in)/1k = 2 mS, so a
    1 V phasor on `in` drives 2 mA from p to n and the 1 kohm loads see
    V(p) = -2 V, V(n) = +2 V.
    Tolerance: gmin (1e-12 S) sits in parallel with the 1 kohm loads, so it
    leaks gmin*V = 2e-12 A -> 2e-9 V across them; measured deviations are
    1.0e-9 V on DC and 2.0e-9 V on AC, so abs=1e-6 (DC) and rel=1e-6 (AC, a
    2e-6 V window) keep 1000x margin and only a wrong stamp can break them.
    """
    results = _run(
        "Behavioral floating I\n"
        "V1 in 0 dc 1 ac 1 0\n"
        "B2 p n I={V(in)**2/1k}\n"
        "R1 p 0 1k\n"
        "R2 n 0 1k\n"
        ".op\n"
        ".ac lin 1 1000 1000\n"
        ".end\n"
    )
    op, ac = results[0], results[1]
    assert op.nodes["p"] == pytest.approx(-1.0, abs=1e-6)
    assert op.nodes["n"] == pytest.approx(1.0, abs=1e-6)
    assert ac.nodes["p"][0] == pytest.approx(-2.0 + 0.0j, rel=1e-6, abs=1e-9)
    assert ac.nodes["n"][0] == pytest.approx(2.0 + 0.0j, rel=1e-6, abs=1e-9)
