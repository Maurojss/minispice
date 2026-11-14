"""Tests for the small-signal AC analysis (.ac)."""

import json
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
