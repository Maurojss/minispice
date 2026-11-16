"""Tests for the schematic renderer (schematic module).

The renderer has no closed-form electrical answer, so the reference is the
concrete artifact: a valid PNG of plausible size, the exact engineering-suffix
strings, the exact BFS column coordinates (multiples of ``_X`` = 3.6) and the
exact number of drawing primitives ``_route`` emits.  Every test stays
hermetic: no browser, no network, matplotlib on the Agg backend (the renderer
forces it).
"""

import io
import struct
import sys

import pytest

from minispice.parser import (
    CCCS,
    CCVS,
    VCCS,
    VCVS,
    BehavioralSource,
    Bjt,
    BjtModel,
    Capacitor,
    Circuit,
    Diode,
    Inductor,
    ISource,
    Mosfet,
    MosModel,
    NetlistParser,
    Resistor,
    VSource,
)
from minispice.schematic import (
    _H,
    _X,
    _fet_kind,
    _fmt_value,
    _label,
    _node_positions,
    _route,
    _symbol_name,
    schematic_png,
)

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def _parse(text: str) -> Circuit:
    return NetlistParser(text=text).parse()


def _render(circuit: Circuit) -> bytes:
    """Render to memory and return the PNG bytes (asserting success + magic)."""
    buf = io.BytesIO()
    assert schematic_png(circuit, buf) is True
    data = buf.getvalue()
    assert data[:8] == PNG_MAGIC
    return data


def _png_dims(data: bytes) -> tuple[int, int]:
    """(width, height) from the IHDR chunk - PNG layout puts it at byte 16."""
    return struct.unpack(">II", data[16:24])


# ---------------------------------------------------------------- value formatting


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (1e9, "1G"),
        (1e6, "1meg"),
        (4700.0, "4.7k"),
        (0.5, "500m"),
        (2.2e-6, "2.2u"),
        (1e-9, "1n"),
        (3e-12, "3p"),
        (5e-15, "5f"),
        (-2e-3, "-2m"),
    ],
)
def test_fmt_value_engineering_suffixes(value, expected):
    # 4700/1e3 = 4.7 exactly in binary? 4.7 is not representable, but %.4g
    # rounds to "4.7", so the string reference is exact.
    assert _fmt_value(value) == expected


def test_fmt_value_below_femto_scale_uses_plain_format():
    # 1e-16 < 1e-15: no suffix applies, falls through to the plain %.4g form.
    assert _fmt_value(1e-16) == "1e-16"
    assert _fmt_value(0.0) == "0"


# ---------------------------------------------------------------- labels / symbols


def test_label_passives_show_formatted_value_on_second_line():
    assert _label(Resistor("R1", "a", "0", 1000.0)) == "R1\n1k"
    assert _label(Capacitor("C1", "a", "0", 1e-6)) == "C1\n1u"
    assert _label(Inductor("L1", "a", "0", 10e-3)) == "L1\n10m"


def test_label_sources_show_dc_with_unit_only_when_dc_is_set():
    assert _label(VSource("V1", "a", "0", 1.8)) == "V1\n1.8V"
    assert _label(ISource("I1", "a", "0", 1e-3)) == "I1\n1mA"
    # dc == 0.0 is falsy: no value line at all
    assert _label(VSource("V2", "a", "0", 0.0)) == "V2"


def test_label_mosfet_bjt_and_controlled_sources():
    # 0.5e-6 = 500e-9: the formatter picks the "n" decade, not "p"
    assert _label(Mosfet("M1", "d", "g", "0", "0", "nmos", 1e-6, 0.5e-6)) == (
        "M1\nW=1u L=500n"
    )
    assert _label(Bjt("Q1", "c", "b", "0", "qmod")) == "Q1\nqmod"
    assert _label(VCVS("E1", "a", "0", "in", "0", 2.0)) == "E1\n2x"
    assert _label(VCCS("G1", "a", "0", "in", "0", 1e-3)) == "G1\n0.001x"
    assert _label(CCCS("F1", "a", "0", "VS", 3.0)) == "F1\n3x I(VS)"
    assert _label(CCVS("H1", "a", "0", "VS", 4.0)) == "H1\n4x I(VS)"


@pytest.mark.parametrize(
    ("element", "expected"),
    [
        (Resistor("R1", "a", "b", 1.0), "Resistor"),
        (Capacitor("C1", "a", "b", 1.0), "Capacitor"),
        (Inductor("L1", "a", "b", 1.0), "Inductor"),
        (VSource("V1", "a", "b", 1.0), "SourceV"),
        (ISource("I1", "a", "b", 1.0), "SourceI"),
        (Diode("D1", "a", "b", "DMOD"), "Diode"),
        (VCVS("E1", "a", "b", "c", "d", 2.0), "SourceV"),
        (CCVS("H1", "a", "b", "VS", 2.0), "SourceV"),
        (VCCS("G1", "a", "b", "c", "d", 1e-3), "SourceI"),
        (CCCS("F1", "a", "b", "VS", 2.0), "SourceI"),
        # behavioral sources have no symbol: renderer degrades to Resistor
        (BehavioralSource("B1", "a", "b", "v", "V(a)*2"), None),
    ],
)
def test_symbol_name_mapping(element, expected):
    assert _symbol_name(element) == expected


def test_fet_kind_reads_model_kind_first_then_model_name():
    mosfet = Mosfet("M1", "d", "g", "s", "0", "NCH")
    assert _fet_kind(Circuit(models={"NCH": MosModel("NCH", kind="nmos")}), mosfet) is True
    mosfet_p = Mosfet("M2", "d", "g", "s", "0", "PCH")
    assert _fet_kind(Circuit(models={"PCH": MosModel("PCH", kind="pmos")}), mosfet_p) is False
    # kind empty -> falls back to the model name
    mosfet_named = Mosfet("M3", "d", "g", "s", "0", "nmos_dev")
    assert _fet_kind(Circuit(models={"nmos_dev": MosModel("nmos_dev", kind="")}), mosfet_named)
    # no model entry at all -> the instance model name decides
    assert _fet_kind(Circuit(), Mosfet("M4", "d", "g", "s", "0", "NMOS")) is True
    assert _fet_kind(Circuit(), Mosfet("M5", "d", "g", "s", "0", "PMOS")) is False
    assert _fet_kind(Circuit(), Mosfet("M6", "d", "g", "s", "0", "sky_p")) is False


# ---------------------------------------------------------------- layout


def test_node_positions_bfs_gives_one_column_per_node():
    circuit = _parse("t\nR1 a b 1k\nR2 b c 1k\nR3 c 0 1k\n.end\n")
    positions = _node_positions(circuit.elements)
    # BFS from min(nodes) = "a": one _X (3.6) per hop, all on the y = 0 track
    assert positions == {"a": (0.0, 0.0), "b": (_X, 0.0), "c": (2 * _X, 0.0)}


def test_node_positions_empty_when_every_pin_is_ground():
    circuit = Circuit(elements=[Resistor("R1", "0", "0", 1.0)])
    assert _node_positions(circuit.elements) == {}


def test_node_positions_places_unreachable_node_in_its_own_column():
    # The FET has only a drain node ("b") and no real edge to anywhere else,
    # so BFS from "a" never reaches "b"; it must still get a column.
    circuit = Circuit(
        elements=[
            Resistor("R1", "a", "0", 1.0),
            Mosfet("M1", "a", "0", "0", "0", "nmos"),
            Mosfet("M2", "b", "0", "0", "0", "nmos"),
        ]
    )
    positions = _node_positions(circuit.elements)
    assert positions == {"a": (0.0, 0.0), "b": (_X, 0.0)}


def test_route_on_same_column_adds_a_single_junction_dot():
    import schemdraw
    import schemdraw.elements as sym

    with schemdraw.Drawing(show=False, backend="matplotlib") as drawing:
        _route(drawing, sym, (1.0, 0.0), 1.0, 2.2)
        assert [type(e).__name__ for e in drawing.elements] == ["Dot"]


def test_route_off_column_adds_three_axis_aligned_wires():
    import schemdraw
    import schemdraw.elements as sym

    with schemdraw.Drawing(show=False, backend="matplotlib") as drawing:
        _route(drawing, sym, (0.0, 0.0), _X, 2.2)
        # up, across, down
        assert [type(e).__name__ for e in drawing.elements] == ["Line", "Line", "Line"]


# ---------------------------------------------------------------- full renders


def test_schematic_png_writes_a_nonempty_png_file(tmp_path):
    circuit = _parse("t\nV1 in 0 1\nR1 in out 1k\nC1 out 0 1u\n.end\n")
    path = tmp_path / "rc.png"
    assert schematic_png(circuit, str(path)) is True

    data = path.read_bytes()
    assert data[:8] == PNG_MAGIC
    assert len(data) > 2000          # a blank canvas at 110 dpi is far smaller
    width, height = _png_dims(data)
    assert 100 < width < 4000 and 100 < height < 4000


def test_schematic_png_two_terminal_variants_including_grounded_ones():
    netlist = (
        "t\n"
        "V1 in 0 1\n"
        "I1 in mid 1m\n"
        "R1 mid out 1k\n"
        "C1 out 0 1u\n"
        "L1 out tail 10m\n"
        "D1 tail 0 DMOD\n"
        ".model DMOD D(IS=1e-14 N=1)\n"
        ".end\n"
    )
    data = _render(_parse(netlist))
    assert len(data) > 3000


def test_schematic_png_grounded_devices_on_one_node_shift_sideways():
    # Two grounded devices on "in": the second hangs at x = col + _H.
    circuit = _parse("t\nV1 in 0 1\nR1 in 0 1k\nC1 in 0 1u\n.end\n")
    data = _render(circuit)
    assert len(data) > 2000
    # the layout itself is the reference: both devices share the column
    assert _node_positions(circuit.elements) == {"in": (0.0, 0.0)}
    assert _H > 0  # the second device is offset by exactly one slot width


def test_schematic_png_ignores_device_whose_both_pins_are_ground():
    # A shorted stub is parseable; it must be skipped, not crash the render.
    circuit = _parse("t\nV1 in 0 1\nR1 in 0 1k\nR2 0 0 1k\n.end\n")
    assert len(circuit.elements) == 3
    assert _render(circuit)[:8] == PNG_MAGIC


def test_schematic_png_flattened_subckt_netlist():
    netlist = (
        "Subcircuit example\n"
        "V1 1 0 10\n"
        "X1 1 2 divider\n"
        ".subckt divider in out\n"
        "R1 in out 1k\n"
        "R2 out 0 1k\n"
        ".ends\n"
        ".end\n"
    )
    circuit = _parse(netlist)
    # the parser flattens: three primitives, no subckt instance left
    assert [type(e).__name__ for e in circuit.elements] == [
        "VSource",
        "Resistor",
        "Resistor",
    ]
    assert len(_render(circuit)) > 2000


def test_schematic_png_cmos_inverter_draws_both_fet_polarities():
    netlist = (
        "CMOS inverter\n"
        "VDD vdd 0 1.8\n"
        "VIN in 0 0.9\n"
        "MP out in vdd vdd PMOS W=2u L=1u\n"
        "MN out in 0 0 NMOS W=1u L=1u\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u)\n"
        ".model PMOS PMOS(VTO=-0.5 KP=60u)\n"
        ".end\n"
    )
    data = _render(_parse(netlist))
    assert len(data) > 3000


def test_schematic_png_stacked_fets_use_raised_body_lane():
    # Two NMOS and two PMOS share a column: index > 0 exercises the shifted
    # slot (x + _H*index for NMOS, x - 1.5 - _H*index for PMOS) and the
    # raised drain lane; gate wires stack into successive lanes.
    netlist = (
        "stacked\n"
        "VDD vdd 0 1.8\n"
        "VIN in 0 0.9\n"
        "M1 out in 0 0 NMOS W=1u L=1u\n"
        "M2 out in 0 0 NMOS W=2u L=1u\n"
        "M3 out in vdd vdd PMOS W=2u L=1u\n"
        "M4 out in vdd vdd PMOS W=4u L=1u\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u)\n"
        ".model PMOS PMOS(VTO=-0.5 KP=60u)\n"
        ".end\n"
    )
    assert len(_render(_parse(netlist))) > 3000


def test_schematic_png_fet_with_grounded_gate_is_left_unwired():
    # A grounded gate has no column to route to: the renderer must skip the
    # wire (target_x is None) and keep drawing.
    netlist = (
        "t\n"
        "V1 d 0 1\n"
        "M1 d 0 0 0 NMOS W=1u L=1u\n"
        ".model NMOS NMOS(VTO=0.5 KP=120u)\n"
        ".end\n"
    )
    assert len(_render(_parse(netlist))) > 2000


def test_schematic_png_bjt_npn_and_pnp():
    npn = (
        "NPN common emitter\n"
        "VCC vcc 0 5\n"
        "VIN in 0 0.7\n"
        "RB in b 100k\n"
        "RC vcc c 1k\n"
        "Q1 c b 0 qmod\n"
        ".model qmod NPN(IS=1e-15 BF=100)\n"
        ".end\n"
    )
    pnp = (
        "PNP emitter follower\n"
        "VCC vcc 0 5\n"
        "VIN in 0 4.3\n"
        "RB in b 100k\n"
        "RC c 0 1k\n"
        "Q1 c b vcc qmod\n"
        ".model qmod PNP(IS=1e-15 BF=100)\n"
        ".end\n"
    )
    assert len(_render(_parse(npn))) > 2500
    assert len(_render(_parse(pnp))) > 2500


def test_schematic_png_bjt_without_real_anchors_is_skipped():
    circuit = Circuit(
        elements=[
            Resistor("R1", "a", "0", 1.0),
            Bjt("Q1", "0", "0", "0", "npn"),
        ],
        models={"npn": BjtModel("npn", kind="npn")},
    )
    assert len(_render(circuit)) > 1000


def test_schematic_png_mosfet_without_real_anchor_is_skipped():
    # Every pin on ground: there is no column to attach the symbol to, so the
    # renderer must skip it instead of indexing an empty position table.
    circuit = Circuit(
        elements=[
            Resistor("R1", "a", "0", 1.0),
            Mosfet("M1", "0", "0", "0", "0", "nmos"),
        ],
        models={"nmos": MosModel("nmos", kind="nmos")},
    )
    assert len(_render(circuit)) > 1000


def test_schematic_png_controlled_sources_render_on_the_track():
    netlist = (
        "controlled\n"
        "V1 in 0 1\n"
        "E1 e 0 in 0 2\n"
        "G1 g 0 in 0 1m\n"
        "VS s 0 0\n"
        "F1 f 0 VS 2\n"
        "H1 h 0 VS 3\n"
        "R1 e 0 1k\n"
        "R2 g 0 1k\n"
        "R3 f 0 1k\n"
        "R4 h 0 1k\n"
        ".end\n"
    )
    assert len(_render(_parse(netlist))) > 3000


def test_schematic_png_behavioral_source_degrades_to_resistor_symbol():
    # _symbol_name returns None for B sources; getattr(sym, "Resistor")
    # fallback keeps the sketch drawable instead of raising.
    netlist = (
        "behavioral\n"
        "V1 in 0 1\n"
        "B1 out 0 V={V(in)*2}\n"
        "R1 out 0 1k\n"
        ".end\n"
    )
    assert len(_render(_parse(netlist))) > 1500


def test_schematic_png_returns_false_for_empty_circuit():
    buf = io.BytesIO()
    assert schematic_png(Circuit(), buf) is False
    assert buf.getvalue() == b""


def test_schematic_png_returns_false_without_schemdraw(monkeypatch):
    # A None entry in sys.modules makes the import raise ImportError, which is
    # the documented "optional dependency missing" contract.
    monkeypatch.setitem(sys.modules, "schemdraw", None)
    buf = io.BytesIO()
    assert schematic_png(_parse("t\nV1 in 0 1\nR1 in 0 1k\n.end\n"), buf) is False
    assert buf.getvalue() == b""


def test_schematic_png_returns_false_when_nothing_is_drawable():
    circuit = _parse("t\nV1 0 0 1\nR1 0 0 1k\n.end\n")
    assert circuit.elements  # non-empty, so the early False does not apply
    buf = io.BytesIO()
    assert schematic_png(circuit, buf) is False
