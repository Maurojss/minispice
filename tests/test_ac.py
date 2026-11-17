"""Tests for the small-signal AC analysis (.ac)."""

import json
import math
import threading
import urllib.request

import numpy as np
import pytest

from minispice import NetlistError, NetlistParser
from minispice.results import AcResult
from minispice.solver import Simulation

LOWPASS = """AC lowpass
V1 in 0 dc 0 ac 1 0
R1 in out 1k
C1 out 0 1u
.ac dec 20 1 100K
.end
"""


# Frequency where omega*R*C = 1000 rad/s * 1 kohm * 1 uF = 1 exactly, so a
# capacitor becomes the pure susceptance j*1 mS and the arithmetic stays
# hand-checkable.
_OMEGA_RC = 1000.0
_F_RC = _OMEGA_RC / (2.0 * math.pi)


def _lowpass():
    circuit = NetlistParser(text=LOWPASS).parse()
    results = Simulation(circuit).run()
    assert len(results) == 1 and isinstance(results[0], AcResult)
    return results[0]


def test_ac_matches_analytic_transfer_function():
    """Every grid point must match H(s) = 1 / (1 + s/fc) to high precision."""
    ac = _lowpass()
    fc = 1.0 / (2 * np.pi * 1e3 * 1e-6)
    imag = 2j * np.pi * ac.freqs
    expected = 1.0 / (1.0 + imag / (2 * np.pi * fc))
    assert np.allclose(ac.nodes["out"], expected, rtol=1e-8)
    assert np.allclose(ac.nodes["in"], np.ones_like(ac.freqs), rtol=1e-9)


def test_ac_bode_values_at_endpoints():
    ac = _lowpass()
    out = ac.nodes["out"]
    db = 20 * np.log10(np.abs(out))
    assert abs(db[0]) < 1e-2                    # passband: 0 dB at fstart
    rolloff = db[-1] + 20 * np.log10(ac.freqs[-1] / (1 / (2 * np.pi * 1e3 * 1e-6)))
    assert abs(rolloff) < 1e-3                  # -20 dB/decade above fc
    # phase reaches ~-90 deg in the stop band and ~0 deg in the passband
    assert abs(np.rad2deg(np.angle(out[0]))) < 1.0
    assert abs(np.rad2deg(np.angle(out[-1])) + 90.0) < 1.0


def test_ac_triggers_at_the_cutoff_frequency():
    """Linear sweep with a grid point exactly at fc."""
    fc = 1.0 / (2 * np.pi * 1e3 * 1e-6)
    netlist = LOWPASS.replace(".ac dec 20 1 100K", f".ac lin 3 {fc:.9g} {2 * fc:.9g}")
    ac = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    assert len(ac.freqs) == 3
    assert ac.freqs[0] == pytest.approx(fc)
    assert abs(abs(ac.nodes["out"][0]) - 1 / np.sqrt(2)) < 1e-6       # -3 dB gain
    assert abs(np.rad2deg(np.angle(ac.nodes["out"][0])) + 45.0) < 1e-4
    db = 20 * np.log10(np.abs(ac.nodes["out"][0]))
    assert abs(db - (-3.0103)) < 1e-3


def test_ac_frequency_kinds_and_counts():
    circuit_text = LOWPASS.rsplit(".ac", 1)[0]
    for kind, points, expected_count in (
        ("dec", 5, 5 * 5 + 1),        # 5 decades of 1..100k -> 25 intervals
        ("oct", 2, 34),               # round(2 * log2(1e5)) = 33 intervals
        ("lin", 11, 11),
    ):
        netlist = f"{circuit_text}.ac {kind} {points} 1 100k\n.end\n"
        ac = Simulation(NetlistParser(text=netlist).parse()).run()[0]
        assert len(ac.freqs) == expected_count
        assert ac.freqs[0] == 1.0 and ac.freqs[-1] == 100e3
        if kind != "lin":
            assert np.all(np.diff(ac.freqs) > 0)


def test_ac_current_source():
    netlist = """tester
I1 0 out dc 0 ac 1m 0
R1 out 0 1k
.ac lin 1 100 100
.end
"""
    ac = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    # 1 mA across 1 kOhm -> 1 V = 0 dB (gmin leaks ~1e-9)
    assert abs(ac.nodes["out"][0] - 1.0) < 1e-8


def test_ac_sources_without_ac_spec_are_short_circuits():
    """A source with no `ac` spec contributes nothing to the AC solution."""
    netlist = """tester
V1 in 0 dc 2
R1 in mid 1k
R2 mid 0 1k
.ac lin 1 1000 1000
.end
"""
    ac = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    assert abs(ac.nodes["in"][0]) < 1e-12
    assert abs(ac.nodes["mid"][0]) < 1e-12


def test_ac_resistive_divider():
    netlist = """divider
V1 in 0 dc 2 ac 2 0
R1 in mid 1k
R2 mid 0 1k
.ac lin 1 1000 1000
.end
"""
    ac = Simulation(NetlistParser(text=netlist).parse()).run()[0]
    assert abs(ac.nodes["mid"][0] - 1.0) < 1e-8          # 2 V / 2 = 1 V
    assert abs(abs(ac.currents["V1"][0]) - 1e-3) < 1e-8  # |2 V - 1 V| / 1 k


def test_ac_inductor_series_rl_complex_phasors():
    """L1 between two non-ground nodes: R1 - L1 - R2 series loop to ground.

    At omega = 1e6 rad/s the 1 mH inductor is j*1 kohm, so the series loop has
    I = 1 / (2k + j1k) = (2 - j)/5 mA = 0.4 - 0.2j mA (Vin = 1 V phasor).
    Then V(a) = 1 - I*1k = 0.6 + 0.2j V, V(b) = I*1k = 0.4 - 0.2j V, and the
    inductor branch current (n1 -> n2) is the same loop current.
    Tolerance: one complex linear solve, so the only error sources are the gmin
    shunt (1e-12 S in parallel with the 1 kohm arms -> a 1e-9 relative offset)
    and float roundoff. Measured deviations are 4.5e-10 V on the 0.6 V phasor
    (7.5e-10 relative) and 2.0e-13 A on the branch current, so rel=1e-6 keeps
    more than 1000x margin: only a wrong stamp can break the test.
    """
    f = 1e6 / (2.0 * math.pi)                    # omega = 1e6 rad/s -> wL = j1 kohm
    ac = Simulation(
        NetlistParser(
            "RL series\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in a 1k\n"
            "L1 a b 1m\n"
            "R2 b 0 1k\n"
            f".ac lin 1 {f:.12f} {f:.12f}\n"
            ".end\n"
        ).parse()
    ).run()[0]
    current = complex(0.4e-3, -0.2e-3)           # (2 - j)/5 mA
    assert ac.nodes["a"][0] == pytest.approx(1.0 - current * 1e3, rel=1e-6, abs=1e-9)
    assert ac.nodes["b"][0] == pytest.approx(current * 1e3, rel=1e-6, abs=1e-9)
    assert ac.nodes["a"][0] == pytest.approx(complex(0.6, 0.2), rel=1e-6, abs=1e-9)
    assert ac.nodes["b"][0] == pytest.approx(complex(0.4, -0.2), rel=1e-6, abs=1e-9)
    assert ac.currents["L1"][0] == pytest.approx(current, rel=1e-6, abs=1e-12)
    # KVL across L1: V(a) - V(b) = j*omega*L*I = j*1000*(0.4-0.2j) mA = 0.2+0.4j V.
    assert ac.nodes["a"][0] - ac.nodes["b"][0] == pytest.approx(
        complex(0.2, 0.4), rel=1e-6, abs=1e-9
    )
    # KCL at `in`: only V1 and R1 meet there -> I(V1) = (V(a) - 1)/1k = -I.
    assert ac.currents["V1"][0] == pytest.approx(-current, rel=1e-6, abs=1e-12)


def test_ac_inductor_rl_divider_transfer_function():
    """RL divider: H(jw) = jwL/(R + jwL) with the inductor leg to ground.

    With R = 1 kohm and L = 1 mH the corner is at omega = R/L = 1e6 rad/s
    (f = 159.154943 kHz). At that frequency wL = R = 1 kohm, so
    H = j/(1 + j) = (1 + j)/2 = 0.5 + 0.5j: magnitude 1/sqrt(2) = -3.0103 dB
    and phase +45 deg, exactly like an RC lowpass at its cutoff. The branch
    current is I = V(out)/(jwL) = (0.5 - 0.5j) mA.

    Tolerance: one complex linear solve; the gmin shunt (1e-12 S against the
    1 mS of R) leaves a 1e-9 relative offset on V(out) (5e-10 V), which is
    5e-13 A in the branch currents; the measured deviation is 1.1e-12 A, so
    abs=1e-11 keeps 10x margin there and rel=1e-9 keeps 1000x on the phasor.
    """
    f = 1e6 / (2.0 * math.pi)                    # omega = R/L = 1e6 rad/s
    ac = Simulation(
        NetlistParser(
            "RL divider\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in out 1k\n"
            "L1 out 0 1m\n"
            f".ac lin 1 {f:.12f} {f:.12f}\n"
            ".end\n"
        ).parse()
    ).run()[0]
    h = complex(0.5, 0.5)                       # jwL/(R + jwL)
    assert ac.nodes["out"][0] == pytest.approx(h, rel=1e-9, abs=1e-12)
    assert abs(ac.nodes["out"][0]) == pytest.approx(1 / math.sqrt(2), rel=1e-9)
    assert math.degrees(math.atan2(ac.nodes["out"][0].imag, ac.nodes["out"][0].real)) \
        == pytest.approx(45.0, abs=1e-6)
    assert 20 * math.log10(abs(ac.nodes["out"][0])) == pytest.approx(-3.0103, abs=1e-3)
    # The inductor leg carries I = V(out)/(jwL) = (0.5 - 0.5j) mA.
    assert ac.currents["L1"][0] == pytest.approx(
        complex(0.5e-3, -0.5e-3), rel=1e-9, abs=1e-11
    )
    # KCL at `out`: the same current leaves through R1, so I(V1) = -0.5+0.5j mA.
    assert ac.currents["V1"][0] == pytest.approx(
        complex(-0.5e-3, 0.5e-3), rel=1e-9, abs=1e-11
    )


def test_ac_voltage_source_between_two_non_ground_nodes():
    """V2 hangs between nodes a and b; neither terminal is node 0.

    V1 drives `in` with a 1 V phasor, R1 in-a, V2 a-b (2 V phasor, 0 V DC),
    R2 b-0. KVL: V(a) - V(b) = 2. Summing the KCL rows at a and b cancels the
    V2 branch current: (V(a) - 1)/1k + V(b)/1k = 0 -> V(a) + V(b) = 1.
    Hence V(a) = 1.5 V, V(b) = -0.5 V, and I(V2) = V(b)/1k = -0.5 mA (the
    branch current is defined a -> b, and only R2 meets node b).
    Tolerance: the only error source is the gmin shunt (1e-12 S against the
    1 mS arms -> 1e-9 relative) plus float roundoff; measured deviations are
    5.0e-10 V on the node voltages and 1.0e-12 A on the branch currents, so
    rel=1e-6 on a 1.5 V phasor (a 1.5e-6 V window) keeps 3000x margin and only
    a wrong stamp can break the test.
    """
    ac = Simulation(
        NetlistParser(
            "Floating V source AC\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in a 1k\n"
            "V2 a b dc 0 ac 2 0\n"
            "R2 b 0 1k\n"
            ".ac lin 1 1000 1000\n"
            ".end\n"
        ).parse()
    ).run()[0]
    assert ac.nodes["a"][0] == pytest.approx(1.5 + 0.0j, rel=1e-6, abs=1e-9)
    assert ac.nodes["b"][0] == pytest.approx(-0.5 + 0.0j, rel=1e-6, abs=1e-9)
    assert ac.nodes["a"][0] - ac.nodes["b"][0] == pytest.approx(
        2.0 + 0.0j, rel=1e-6, abs=1e-9
    )
    assert ac.currents["V2"][0] == pytest.approx(-0.5e-3 + 0.0j, rel=1e-6, abs=1e-12)
    assert ac.currents["V1"][0] == pytest.approx(0.5e-3 + 0.0j, rel=1e-6, abs=1e-12)


def test_ac_current_source_into_rc_load_from_non_ground_node():
    """I1 pulls 0.5 mA out of the non-ground node `a`; `a` also sees R and C.

    With R1 = 1 kohm from the 1 V source, R2 = 1 kohm and C1 = 1 uF to ground,
    at omega*R*C = 1 the parallel load is 1 mS + j*1 mS. I1 sinks 0.5 mA out of
    `a`, so KCL gives V(a)*(2 + j) mS = 1 mA - 0.5 mA and
    V(a) = 0.5/(2 + j) = 0.2 - 0.1j V (the I/Z closed form).
    Tolerance: gmin adds 1e-12 S to the 2 mS sum (~5e-10 relative), so
    rel=1e-6 / abs=1e-9 catches a sign flip in the b[i] stamp without
    flakiness (measured deviation 1.0e-10 V on the phasor).
    """
    ac = Simulation(
        NetlistParser(
            "Current source AC\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in a 1k\n"
            "R2 a 0 1k\n"
            "C1 a 0 1u\n"
            "I1 a 0 dc 0 ac 0.5m 0\n"
            f".ac lin 1 {_F_RC:.12f} {_F_RC:.12f}\n"
            ".end\n"
        ).parse()
    ).run()[0]
    assert ac.nodes["a"][0] == pytest.approx(complex(0.2, -0.1), rel=1e-6, abs=1e-9)
    # KCL at `in`: only V1 and R1 meet there, so I(V1) = (V(a) - 1)/1k.
    assert ac.currents["V1"][0] == pytest.approx(
        (complex(0.2, -0.1) - 1.0) / 1e3, rel=1e-6, abs=1e-12
    )


@pytest.mark.parametrize("temp_c", [27.0, 100.0])
def test_ac_diode_small_signal_uses_operating_point_conductance(temp_c):
    """A DC-biased diode is a conductance gd = Id/rd with rd = n*Vt/Id in AC.

    I1 forces 1 mA through D1 (C1 blocks DC, so the diode carries the whole
    bias), which sets Vd = n*Vt*ln(1 + Id/IS) and the small-signal resistance
    rd = n*Vt/Id. The .op report carries `gd` (Id/rd) and `id`; from the
    reported gd the AC divider of C1 = 1 nF at 1 MHz against the diode is the
    closed form V(a) = j*wC/(gd + j*wC).

    The model pins TNOM to the netlist temperature, so the IS(T) scaling
    (T/TNOM)^(XTI/N)*exp(EG/(N*Vt)*(T/TNOM-1)) collapses to 1 and the
    saturation current used by the solver is exactly the model's IS. Vt =
    kT/q is evaluated at the .temp of the netlist, which is what makes the
    reference hand-derivable at both temperatures.

    Tolerance: the AC solve is a single complex linear system, so the
    reference built from the reported gd matches to gmin/gd = 2.6e-11 at 27 C
    and 3.2e-11 at 100 C (measured 2.6e-11 and 3.2e-11); rel=1e-9 keeps 40x
    margin. The DC point itself is only accurate to Newton's stopping rule
    (|dV| <= vntol + reltol*|V| = 1e-6 + 1e-3*0.655 = 656 uV), so Vd is
    checked with abs=1e-3 and gd and rd with rel=5e-2 (measured 4.6e-6 V and
    1.8e-4 at 27 C, 6.5e-8 V and 2.0e-6 at 100 C, i.e. at least 140x inside).
    """
    vt = 1.380649e-23 * (temp_c + 273.15) / 1.602176634e-19   # kT/q
    idc = 1e-3                                               # forced bias
    isat = 1e-14                                             # model IS
    n = 1.0                                                  # model N
    vd_analytic = n * vt * math.log(1 + idc / isat)
    rd_analytic = n * vt / idc                               # 25.85 ohm at 27 C
    gd_analytic = 1.0 / rd_analytic                           # Id/rd = 38.67 mS
    omega_c = 2.0 * math.pi * 1e6 * 1e-9                    # wC = 6.283 mS
    results = Simulation(
        NetlistParser(
            "Diode AC\n"
            "V1 s 0 dc 0 ac 1 0\n"
            "C1 s a 1n\n"
            "I1 0 a dc 1m\n"
            "D1 a 0 dmod\n"
            f".model dmod D(IS=1e-14 N=1 TNOM={temp_c})\n"
            f".temp {temp_c}\n"
            ".op\n"
            ".ac lin 1 1e6 1e6\n"
            ".end\n"
        ).parse()
    ).run()
    op, ac = results[0], results[1]
    assert op.nodes["a"] == pytest.approx(vd_analytic, abs=1e-3)
    assert op.details["D1"]["id"] == pytest.approx(idc, rel=1e-2)
    assert op.details["D1"]["gd"] == pytest.approx(gd_analytic, rel=5e-2)
    # rd = n*Vt/Id must be the reciprocal of the reported conductance.
    assert 1.0 / op.details["D1"]["gd"] == pytest.approx(rd_analytic, rel=5e-2)
    gd = op.details["D1"]["gd"]
    assert ac.nodes["a"][0] == pytest.approx(
        1j * omega_c / (gd + 1j * omega_c), rel=1e-9, abs=1e-12
    )
    # KCL at `s`: the only path is C1, so I(V1) = -j*wC*(1 - V(a)).
    assert ac.currents["V1"][0] == pytest.approx(
        -1j * omega_c * (1.0 - ac.nodes["a"][0]), rel=1e-9, abs=1e-12
    )


def test_ac_decade_sweep_with_equal_endpoints_keeps_a_single_point():
    """`.ac dec 10 1000 1000` spans zero decades, so the grid holds one point.

    Closed form: the frequency list is [1000 Hz] and the 1 kohm / 1 kohm
    divider of the 1 V phasor gives V(out) = 0.5 V.
    """
    ac = Simulation(
        NetlistParser(
            "Degenerate sweep\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".ac dec 10 1000 1000\n"
            ".end\n"
        ).parse()
    ).run()[0]
    assert list(ac.freqs) == [1000.0]
    assert ac.nodes["out"][0] == pytest.approx(0.5 + 0.0j, rel=1e-9, abs=1e-12)


def test_ac_zero_resistance_is_rejected():
    """A 0-ohm resistor has no finite AC stamp (1/R would blow up).

    Every AC sweep solves the DC operating point first, and the DC assembly
    raises the same error, so the guard that the netlist triggers is the DC
    one; solver.py keeps a defensive copy of the check in the AC stamping path
    that the DC guard makes unreachable (reported to spice-solver-dev).
    """
    netlist = """zero
V1 in 0 dc 0 ac 1 0
R1 in out 0
R2 out 0 1k
.ac lin 1 1000 1000
.end
"""
    with pytest.raises(ValueError, match="zero resistance"):
        Simulation(NetlistParser(text=netlist).parse()).run()


def test_ac_nonlinear_devices_do_not_crash():
    """Diodes and MOSFETs are linearized around the DC operating point."""
    netlist = """CMOS inverter AC
VDD vdd 0 5
VIN in 0 dc 2.5 ac 1 0
.model n1 nmos kp=1e-4 vt=1 w=2u l=1u
.model p1 pmos kp=1e-4 vt=-1 w=4u l=1u
MN out in 0 0 n1 m=1
MP out in vdd 0 p1 m=1
Cload out 0 1p
.ac dec 5 1 1G
.end
"""
    circuit = NetlistParser(text=netlist).parse()
    ac = Simulation(circuit).run()[0]
    assert np.all(np.isfinite(ac.nodes["out"]))
    assert len(ac.freqs) == int(round(5 * 9)) + 1  # 9 decades, 5/decade


def test_ac_parser_errors():
    base = LOWPASS.rsplit(".ac", 1)[0]
    with pytest.raises(NetlistError, match="dec, oct or lin"):
        NetlistParser(text=base + ".ac log 10 1 100\n.end\n").parse()
    with pytest.raises(NetlistError, match="fstart and fstop"):
        NetlistParser(text=base + ".ac dec 10 1\n.end\n").parse()
    with pytest.raises(NetlistError, match="non-zero"):
        NetlistParser(text=base + ".ac dec 10 0 100\n.end\n").parse()


def test_ac_layer_errors():
    circuit = NetlistParser(text=LOWPASS).parse()
    analysis = circuit.analyses[0]
    simulation = Simulation(circuit)
    with pytest.raises(ValueError, match="fstart"):
        simulation.ac_sweep(analysis.kind, 10, 0.0, 100.0)
    with pytest.raises(ValueError, match="kind"):
        simulation.ac_sweep("bogus", 10, 1.0, 100.0)


def test_ac_result_str_and_plot(tmp_path):
    """The classroom displays and CLI --plot work for AC too."""
    ac = _lowpass()
    text = str(ac)
    assert "AC analysis" in text and "|V(out)|dB" in text
    path = tmp_path / "ac.png"
    ac.plot(nodes=["out"], path=str(path))
    assert path.exists() and path.read_bytes()[:4] == b"\x89PNG"


def test_ac_webui_payload_plot_and_cursor():
    """HTTP flow: /run, /plot (dB, log-x), /dat (log10 x, dB, disp), /csv."""
    server = None
    from minispice.webui import create_server

    server = create_server()
    threading.Thread(target=server.serve_forever, daemon=True).start()
    host, port = server.server_address[:2]
    base = f"http://{host}:{port}"
    try:
        run_req = urllib.request.Request(
            f"{base}/run",
            data=json.dumps({"netlist": LOWPASS}).encode(),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        data = json.loads(urllib.request.urlopen(run_req, timeout=15).read())
        assert data["ok"]
        run = data["run"]
        info = data["results"][0]
        assert info["type"] == "series" and info["title"] == "AC analysis"
        assert info["head"][0] == "freq (Hz)"

        # plot renders a PNG (dB magnitude on a log axis)
        with urllib.request.urlopen(
            f"{base}/plot?run={run}&i=0&sigs=V:out", timeout=15
        ) as response:
            assert response.read(4) == b"\x89PNG"
        # zoom works in Hz
        with urllib.request.urlopen(
            f"{base}/plot?run={run}&i=0&sigs=V:out&x0=50&x1=2e4", timeout=15
        ) as response:
            assert response.read(4) == b"\x89PNG"

        # /dat: values already in dB, x in log10(f), disp with suffix
        with urllib.request.urlopen(f"{base}/dat?run={run}&i=0", timeout=15) as response:
            dat = json.loads(response.read())
        assert dat["axis"] == "freq (Hz)"
        assert len(dat["disp"]) == len(dat["x"]) == len(dat["series"]["V:out"])
        assert dat["x"][0] == 0.0                       # log10(1 Hz)
        assert dat["disp"][-1] == "100k" + "Hz"
        fc = 1 / (2 * np.pi * 1e3 * 1e-6)
        nearest = int(np.argmin(np.abs(np.asarray(dat["x"]) - np.log10(fc))))
        assert abs(dat["series"]["V:in"][nearest]) < 0.05       # ~0 dB in-band
        assert -3.6 < dat["series"]["V:out"][nearest] < -2.4     # ~-3 dB at fc

        # /csv exports dB and phase columns
        with urllib.request.urlopen(f"{base}/csv?run={run}&i=0", timeout=15) as response:
            csv_text = response.read().decode()
        lines = csv_text.splitlines()
        assert lines[0] == "freq,V:in dB,V:out dB,I:V1 dB,V:in deg,V:out deg,I:V1 deg"
        assert lines[-1].split(",")[0].startswith("100000.0")
    finally:
        server.shutdown()
        server.server_close()
