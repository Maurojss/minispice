import pytest

from minispice import NetlistParser
from minispice.parser import NetlistError, PulseWave, SinWave, parse_value


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("10", 10.0),
        ("1k", 1e3),
        ("2.2u", 2.2e-6),
        ("1meg", 1e6),
        ("1M", 1e-3),
        ("4.7n", 4.7e-9),
        ("1e-3", 1e-3),
        ("100nF", 100e-9),
        ("10kohm", 10e3),
    ],
)
def test_parse_value(text, expected):
    assert parse_value(text) == pytest.approx(expected)


def test_invalid_value():
    with pytest.raises(NetlistError):
        parse_value("abc")


def test_title_detection():
    circuit = NetlistParser(
        "RC low-pass filter\nV1 in 0 1\nR1 in 0 1k\n.op\n.end\n"
    ).parse()
    assert circuit.title == "RC low-pass filter"
    assert len(circuit.elements) == 2


def test_no_title_netlist():
    circuit = NetlistParser("V1 in 0 1\nR1 in 0 1k\n.op\n.end\n").parse()
    assert circuit.title == ""
    assert len(circuit.elements) == 2


def test_no_title_starting_with_diode():
    circuit = NetlistParser(
        "D1 in 0 dmod\nV1 in 0 1\n.model dmod D(IS=1e-14)\n.op\n.end\n"
    ).parse()
    assert circuit.title == ""
    assert len(circuit.elements) == 2


def test_source_waveforms():
    circuit = NetlistParser(
        "V1 in 0 PULSE(0 1 0 1n 1n 1m 2m)\n"
        "V2 in 0 SIN(0 5 1k)\n"
        "R1 in 0 1k\n.op\n.end\n"
    ).parse()
    pulse = circuit.elements[0]
    sine = circuit.elements[1]
    assert isinstance(pulse.waveform, PulseWave)
    assert isinstance(sine.waveform, SinWave)
    assert pulse.value_at(0.5e-3) == pytest.approx(1.0)
    assert sine.value_at(0.25e-3) == pytest.approx(5.0)


def test_model_params():
    circuit = NetlistParser(
        "M1 d g 0 0 NMOS W=2u L=1u\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)\n"
        ".op\n.end\n"
    ).parse()
    model = circuit.models["nmos"]
    assert model.kind == "nmos"
    assert model.vto == pytest.approx(0.5)
    assert model.kp == pytest.approx(120e-6)
    assert model.lam == pytest.approx(0.02)


def test_subcircuit_flattening():
    circuit = NetlistParser(
        "V1 1 0 10\n"
        "X1 1 2 divider\n"
        ".subckt divider in out\n"
        "R1 in out 1k\n"
        "R2 out 0 1k\n"
        ".ends\n"
        ".op\n.end\n"
    ).parse()
    names = [el.name for el in circuit.elements]
    assert "X1:R1" in names
    assert circuit.node_names == ["1", "2"]


def test_nested_subcircuits():
    circuit = NetlistParser(
        "V1 1 0 10\n"
        "X1 1 2 divider\n"
        ".subckt divider in out\n"
        "X2 in out half\n"
        ".ends\n"
        ".subckt half in out\n"
        "R1 in out 1k\n"
        "R2 out 0 1k\n"
        ".ends\n"
        ".op\n.end\n"
    ).parse()
    names = [el.name for el in circuit.elements]
    assert "X1:X2:R1" in names


def test_include(tmp_path):
    sub = tmp_path / "part.cir"
    sub.write_text("R1 1 0 1k\n")
    main = tmp_path / "main.cir"
    main.write_text("V1 1 0 1\n.include part.cir\n.op\n.end\n")
    circuit = NetlistParser.from_file(main).parse()
    assert len(circuit.elements) == 2
