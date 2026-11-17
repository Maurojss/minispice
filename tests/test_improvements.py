"""Tests for parameters, temperature, initial conditions and measurements."""

import csv
import math

import numpy as np
import pytest

from minispice import NetlistParser, Simulation
from minispice.cli import main as cli_main
from minispice.parser import NetlistError, TranAnalysis, parse_value


def test_param_expressions_and_forward_references():
    circuit = NetlistParser(
        "Param divider\n"
        ".param R2={RV*2}\n"
        ".param RV=1k\n"
        ".param VDD={2*0.5}\n"
        "V1 in 0 {VDD}\n"
        "R1 in out {RV}\n"
        "R2 out 0 {R2}\n"
        ".op\n.end\n"
    ).parse()
    assert circuit.params["rv"] == pytest.approx(1e3)
    assert circuit.params["r2"] == pytest.approx(2e3)
    assert circuit.params["vdd"] == pytest.approx(1.0)
    assert Simulation(circuit).operating_point().nodes["out"] == pytest.approx(2 / 3)


def test_param_inside_model_and_waveform():
    circuit = NetlistParser(
        "Param model\n"
        ".param VTH=0.5\n"
        ".param AMP=2\n"
        "V1 g 0 SIN(0 {AMP} 1k)\n"
        "M1 d g 0 0 NMOS W=2u L=1u\n"
        "VDD d 0 1.8\n"
        ".model NMOS NMOS(VTO={VTH} KP=120u)\n"
        ".op\n.end\n"
    ).parse()
    assert circuit.models["nmos"].vto == pytest.approx(0.5)
    assert circuit.elements[0].waveform.amplitude == pytest.approx(2.0)


def test_invalid_expression_is_reported():
    with pytest.raises(NetlistError, match="unknown parameter"):
        parse_value("{MISSING*2}", {})


def test_temp_directive_changes_the_operating_point():
    netlist = (
        "Temp diode\n"
        "I1 0 out 1m\n"
        "D1 out 0 dmod\n"
        ".model dmod D(IS=1e-14 N=1)\n"
        "{temp}\n"
        ".op\n.end\n"
    )
    cold = NetlistParser(netlist.format(temp=".temp 27")).parse()
    hot = NetlistParser(netlist.format(temp=".temp 100")).parse()
    assert cold.options["temp"] == "27.0"
    assert hot.options["temp"] == "100.0"
    v_cold = Simulation(cold).operating_point().nodes["out"]
    v_hot = Simulation(hot).operating_point().nodes["out"]
    assert v_hot < v_cold


def test_tran_uic_and_capacitor_ic():
    result = Simulation(
        NetlistParser(
            "UIC capacitor\n"
            "R1 out 0 1k\n"
            "C1 out 0 1u IC=0.5\n"
            ".tran 10u 1m uic\n.end\n"
        ).parse()
    ).run()[0]
    assert result.nodes["out"][0] == pytest.approx(0.5)
    assert result.nodes["out"][-1] == pytest.approx(0.5 / math.e, abs=2e-3)


def test_tran_ic_directive_sets_node_voltage():
    circuit = NetlistParser(
        "IC node\n"
        "R1 out 0 1k\n"
        "C1 out 0 1u\n"
        ".ic V(out)=0.5\n"
        ".tran 10u 1m uic\n.end\n"
    ).parse()
    assert isinstance(circuit.analyses[0], TranAnalysis)
    assert circuit.analyses[0].uic is True
    result = Simulation(circuit).run()[0]
    assert result.nodes["out"][0] == pytest.approx(0.5)
    assert result.nodes["out"][-1] == pytest.approx(0.5 / math.e, abs=2e-3)


def test_measurement_delay_and_statistics():
    result = Simulation(
        NetlistParser(
            "Measure RC\n"
            "V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n"
            "R1 in out 1k\n"
            "C1 out 0 1u\n"
            ".tran 10u 5m\n"
            ".meas tran trise TRIG v(out) VAL=0.1 RISE=1 "
            "TARG v(out) VAL=0.9 RISE=1\n"
            ".meas tran vmax MAX v(out)\n"
            ".meas tran vmin MIN v(out)\n"
            ".meas tran vavg AVG v(out)\n"
            ".meas tran vrms RMS v(out)\n"
            ".end\n"
        ).parse()
    ).run()[0]
    assert result.measurements["trise"] == pytest.approx(1e-3 * math.log(9), abs=2e-5)
    assert result.measurements["vmax"] == pytest.approx(1 - math.exp(-5), abs=2e-3)
    assert result.measurements["vmin"] == pytest.approx(0.0, abs=1e-6)
    assert result.measurements["vavg"] == pytest.approx(0.8013, abs=2e-3)
    assert result.measurements["vrms"] == pytest.approx(0.8383, abs=2e-3)


def test_measurement_on_dc_sweep():
    result = Simulation(
        NetlistParser(
            "DC measure\n"
            "V1 in 0 0\n"
            "R1 in out 1k\n"
            "R2 out 0 1k\n"
            ".dc V1 0 1 0.1\n"
            ".meas dc voutmax MAX V(out)\n"
            ".end\n"
        ).parse()
    ).run()[0]
    assert result.measurements["voutmax"] == pytest.approx(0.5)


def test_measurement_of_a_branch_current_signal():
    """`.meas` can measure a branch current, not only a node voltage.

    Closed form: the divider V1 1 0 10 / R1 1 2 1k / R2 2 0 1k draws
    I(V1) = -V1/(R1+R2) = -V1/2k, so over the sweep V1 = 0, 5, 10 the branch
    current is 0, -2.5 mA and -5 mA. MIN is therefore -5 mA and AVG -2.5 mA.
    The bare name `V1` is not a node name (the nodes are `1` and `2`), so the
    lookup falls back to the branch-current map and measures the same value.

    Tolerance: the gmin shunt (1e-12 S) across the divider leaks 1e-12 A at
    the 1 V node, which shows up as a 1.25e-11 A offset (measured); rel=1e-6
    (a 5e-9 A window) keeps 400x margin.

    An unknown signal is reported as NaN plus a warning instead of aborting the
    run, so the remaining measurements of the same card still land.
    """
    circuit = NetlistParser(
        "DC measure current\n"
        "V1 1 0 10\n"
        "R1 1 2 1k\n"
        "R2 2 0 1k\n"
        ".dc V1 0 10 5\n"
        ".meas dc imin MIN I(V1)\n"
        ".meas dc iavg AVG I(V1)\n"
        ".meas dc ibare AVG V1\n"
        ".meas dc ibogus MAX NOPE\n"
        ".end\n"
    ).parse()
    result = Simulation(circuit).run()[0]
    assert list(result.sweep) == pytest.approx([0.0, 5.0, 10.0])
    assert result.measurements["imin"] == pytest.approx(-5e-3, rel=1e-6, abs=1e-9)
    assert result.measurements["iavg"] == pytest.approx(-2.5e-3, rel=1e-6, abs=1e-9)
    assert result.measurements["ibare"] == pytest.approx(-2.5e-3, rel=1e-6, abs=1e-9)
    assert math.isnan(result.measurements["ibogus"])
    assert any("unknown signal in .meas: NOPE" in w for w in circuit.warnings)


def test_operating_point_details_for_mosfets_and_diodes():
    result = Simulation(
        NetlistParser(
            "Details\n"
            "VDD vdd 0 1.8\n"
            "VIN in 0 0.9\n"
            "MP out in vdd vdd PMOS W=2u L=1u\n"
            "MN out in 0 0 NMOS W=1u L=1u\n"
            ".model NMOS NMOS(VTO=0.5 KP=120u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
            ".model PMOS PMOS(VTO=-0.5 KP=60u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
            ".op\n.end\n"
        ).parse()
    ).operating_point()
    assert set(result.details) == {"MP", "MN"}
    assert result.details["MN"]["region"] in {"cutoff", "subthreshold", "linear", "saturation"}
    assert result.details["MN"]["vth"] == pytest.approx(0.5, abs=0.15)
    assert result.details["MP"]["type"] == "pmos"


def test_ac_phase_helper():
    result = Simulation(
        NetlistParser(
            "AC\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "R1 in out 1k\n"
            "C1 out 0 1u\n"
            ".ac dec 20 1 100k\n.end\n"
        ).parse()
    ).run()[0]
    assert result.phase(result.nodes["out"])[0] == pytest.approx(0.0, abs=1.0)
    assert result.phase(result.nodes["out"])[-1] == pytest.approx(-90.0, abs=1.0)


def test_cli_ac_csv_exports_phase(tmp_path):
    netlist = tmp_path / "ac.cir"
    netlist.write_text(
        "AC\nV1 in 0 dc 0 ac 1 0\nR1 in out 1k\nC1 out 0 1u\n"
        ".ac dec 5 10 100k\n.end\n"
    )
    out = tmp_path / "ac.csv"
    assert cli_main([str(netlist), "--csv", str(out), "--quiet"]) == 0
    with out.open() as fh:
        rows = list(csv.reader(fh))
    assert any("deg" in head for head in rows[0])
    assert any("dB" in head for head in rows[0])


def test_bjt_ebers_moll_operating_point():
    result = Simulation(
        NetlistParser(
            "NPN common emitter\n"
            "VCC vcc 0 5\n"
            "VIN in 0 0.7\n"
            "RB in b 100k\n"
            "RC vcc c 1k\n"
            "Q1 c b 0 qmod\n"
            ".model qmod NPN(IS=1e-15 BF=100)\n"
            ".op\n.end\n"
        ).parse()
    ).operating_point()
    ib = (0.7 - result.nodes["b"]) / 100e3
    ic = result.currents["Q1"]
    assert ic == pytest.approx(100 * ib, rel=1e-3)
    assert result.nodes["c"] == pytest.approx(5.0 - ic * 1e3, rel=1e-6)
    assert result.details["Q1"]["region"] == "forward active"
    assert result.details["Q1"]["beta"] == pytest.approx(100.0, rel=1e-2)


def test_dependent_sources_dc_and_ac():
    result = Simulation(
        NetlistParser(
            "Dependent sources\n"
            "V1 in 0 1\n"
            "R0 in 0 1k\n"
            "E1 eout 0 in 0 2\n"
            "R1 eout 0 1k\n"
            "G1 0 gout in 0 1m\n"
            "R2 gout 0 1k\n"
            "F1 fout 0 V1 2\n"
            "R3 fout 0 1k\n"
            "H1 hout 0 V1 2k\n"
            "R4 hout 0 1k\n"
            ".op\n.end\n"
        ).parse()
    ).operating_point()
    assert result.nodes["eout"] == pytest.approx(2.0)
    assert result.nodes["gout"] == pytest.approx(1.0)
    assert result.nodes["fout"] == pytest.approx(2.0)
    assert result.nodes["hout"] == pytest.approx(-2.0)
    assert result.currents["E1"] == pytest.approx(-2e-3)

    ac = Simulation(
        NetlistParser(
            "AC dependent\n"
            "V1 in 0 dc 0 ac 1 0\n"
            "E1 out 0 in 0 3\n"
            "R1 out 0 1k\n"
            ".ac lin 3 10 100\n.end\n"
        ).parse()
    ).run()[0]
    assert np.allclose(ac.nodes["out"], 3.0 + 0j)


def test_nested_dc_is_reported_instead_of_silently_ignored():
    circuit = NetlistParser(
        "Nested\nV1 a 0 0\nV2 b 0 0\nR1 a out 1k\nR2 out b 1k\n"
        ".dc V1 0 1 0.5 V2 0 1 0.5\n.end\n"
    ).parse()
    assert any("nested .dc" in warning for warning in circuit.warnings)
