"""Circuit schematic rendering for the web GUI.

Draws the parsed (already flattened) netlist as a teaching-quality sketch
with schemdraw. Real (non-ground) nodes are laid out on a horizontal track
in BFS order:

* two-terminal devices between two real nodes are drawn straight on the
  track;
* devices with one pin on ground hang vertically from their node column
  down to a ground symbol;
* MOSFETs are drawn horizontally on the track (drain left, source right,
  gate up). Gate wires route through lanes *above* the track; wires that
  would run through another symbol dive into lanes *below* the track.

Every wire is axis-aligned — a readable wiring sketch, not place-and-route.
"""

from __future__ import annotations

import io
from collections import defaultdict, deque

from .parser import (
    CCCS,
    CCVS,
    VCCS,
    VCVS,
    Bjt,
    Capacitor,
    Circuit,
    Diode,
    Inductor,
    ISource,
    Mosfet,
    Resistor,
    VSource,
)

_X = 3.6    # horizontal spacing between node columns
_H = 2.0    # column slot width for devices attached to one node
_HANG = 1.9  # vertical length of a hanging (grounded) device
_GATE_LANE = 2.2     # y of the first gate wire lane
_GATE_STEP = 0.7
_BODY_LANE = 1.15    # y of the raised drain-wire lane (slot > 0)
_WIRE_LANE = -0.9    # y of the first below-track wire lane
_WIRE_STEP = 0.7


def _fmt_value(value: float) -> str:
    scale = [
        (1e9, "G"), (1e6, "meg"), (1e3, "k"), (1.0, ""), (1e-3, "m"),
        (1e-6, "u"), (1e-9, "n"), (1e-12, "p"), (1e-15, "f"),
    ]
    for factor, suffix in scale:
        if abs(value) >= factor:
            return f"{value / factor:.4g}".rstrip(".") + suffix
    return f"{value:.4g}"


def _pin_nodes(element) -> list[tuple[str, str]]:
    """(role, node) pairs (FET bulk joins the source and is not drawn)."""
    if isinstance(element, Mosfet):
        return [("drain", element.nd), ("gate", element.ng), ("source", element.ns)]
    if isinstance(element, Bjt):
        return [
            ("collector", element.nc),
            ("base", element.nb),
            ("emitter", element.ne),
        ]
    return [("n1", element.n1), ("n2", element.n2)]


def _symbol_name(element) -> str | None:
    if isinstance(element, Resistor):
        return "Resistor"
    if isinstance(element, Capacitor):
        return "Capacitor"
    if isinstance(element, Inductor):
        return "Inductor"
    if isinstance(element, VSource):
        return "SourceV"
    if isinstance(element, ISource):
        return "SourceI"
    if isinstance(element, Diode):
        return "Diode"
    if isinstance(element, (VCVS, CCVS)):
        return "SourceV"
    if isinstance(element, (VCCS, CCCS)):
        return "SourceI"
    return None


def _label(element) -> str:
    if isinstance(element, (Resistor, Capacitor, Inductor)):
        value = f"\n{_fmt_value(element.value)}"
    elif isinstance(element, (VSource, ISource)):
        value = f"\n{_fmt_value(element.dc)}{'V' if isinstance(element, VSource) else 'A'}" if element.dc else ""
    elif isinstance(element, Mosfet):
        value = f"\nW={_fmt_value(element.w)} L={_fmt_value(element.l)}"
    elif isinstance(element, Bjt):
        value = f"\n{element.model}"
    elif isinstance(element, (VCVS, VCCS)):
        value = f"\n{element.gain:g}x"
    elif isinstance(element, (CCCS, CCVS)):
        value = f"\n{element.gain:g}x I({element.control})"
    else:
        value = ""
    return f"{element.name}{value}"


def _fet_kind(circuit: Circuit, element: Mosfet) -> bool:
    """True when the device is an NMOS."""
    model = circuit.models.get(element.model)
    kind = getattr(model, "kind", "") if model is not None else ""
    name = (kind or element.model).lower()
    return not (name.startswith("p") or name.endswith("pmos") or "_p" in name)


def _node_positions(elements) -> dict[str, tuple[float, float]]:
    """BFS layout of the real (non-ground) nodes on the y = 0 track."""
    nodes: set[str] = set()
    for element in elements:
        nodes.update(node for _, node in _pin_nodes(element) if node != "0")
    if not nodes:
        return {}

    adjacency: dict[str, set[str]] = {node: set() for node in nodes}
    for element in elements:
        real = [node for _, node in _pin_nodes(element) if node != "0"]
        for a, b in zip(real, real[1:]):
            adjacency[a].add(b)
            adjacency[b].add(a)

    start = min(adjacency)
    positions = {start: (0.0, 0.0)}
    queue = deque([start])
    x = 0.0
    while queue:
        node = queue.popleft()
        for neighbor in sorted(adjacency[node]):
            if neighbor not in positions:
                x += _X
                positions[neighbor] = (x, 0.0)
                queue.append(neighbor)
    for node in sorted(nodes - set(positions)):
        x += _X
        positions[node] = (x, 0.0)
    return positions


def _route(drawing, sym, start, target_x: float, level: float) -> None:
    """Route (start) -> (start, level) -> (target_x, level) -> (target_x, 0)."""
    if abs(start[0] - target_x) < 1e-9:
        drawing.add(sym.Dot().at(start))
        return
    drawing.add(sym.Line().at(start).to((start[0], level)))
    drawing.add(sym.Line().at((start[0], level)).to((target_x, level)))
    drawing.add(sym.Line().at((target_x, level)).to((target_x, 0.0)))


def schematic_png(circuit: Circuit, path: str | io.BytesIO) -> bool:
    """Render the circuit to `path`; True on success, False when an
    optional dependency is missing or nothing is drawable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import schemdraw
        import schemdraw.elements as sym
    except ImportError:
        return False

    elements = list(circuit.elements)
    if not elements:
        return False
    positions = _node_positions(elements)
    hangs: defaultdict[str, int] = defaultdict(int)     # grounded devices per node
    n_fets: defaultdict[str, int] = defaultdict(int)    # NMOS per node column
    p_fets: defaultdict[str, int] = defaultdict(int)    # PMOS per node column
    gate_lanes = 0                                      # gate wires routed above
    wire_lanes = 0                                      # wires routed below
    labeled: set[str] = set()

    with schemdraw.Drawing(show=False, backend="matplotlib") as drawing:
        drawing.config(fontsize=10, lw=1.5)

        def dot(node_name: str, loc: str = "bot") -> None:
            point = positions[node_name]
            drawing.add(sym.Dot().at(point))
            if node_name not in labeled:
                labeled.add(node_name)
                drawing.add(
                    sym.Label().at(point).label(node_name, loc=loc, fontsize=9, color="#555")
                )

        for element in elements:
            label = _label(element)

            if isinstance(element, Bjt):
                model = circuit.models.get(element.model)
                npn = getattr(model, "kind", "npn") != "pnp"
                bjt_cls = sym.BjtNpn if npn else sym.BjtPnp
                anchor_node = next(
                    (
                        node
                        for _, node in _pin_nodes(element)
                        if node != "0" and node in positions
                    ),
                    None,
                )
                if anchor_node is None:
                    continue
                col = positions[anchor_node][0]
                device = drawing.add(
                    bjt_cls().at((col, 0.0)).theta(90).label(label, loc="bot")
                )
                anchors = device.absanchors
                for role, node_name in _pin_nodes(element):
                    anchor_point = anchors.get(role)
                    if node_name == "0":
                        drawing.add(sym.Ground().at(anchor_point))
                        continue
                    target_x = positions.get(node_name, (None,))[0]
                    if target_x is None or anchor_point is None:
                        continue
                    if node_name != anchor_node:
                        level = _WIRE_LANE - _WIRE_STEP * wire_lanes
                        wire_lanes += 1
                        _route(drawing, sym, anchor_point, target_x, level)
                    dot(node_name)
                continue

            if isinstance(element, Mosfet):
                nmos = _fet_kind(circuit, element)
                fet_cls = sym.NFet if nmos else sym.PFet
                anchor_node = next(
                    (node for _, node in _pin_nodes(element) if node != "0" and node in positions),
                    None,
                )
                if anchor_node is None:
                    continue
                col = positions[anchor_node][0]
                index = (n_fets if nmos else p_fets)[anchor_node]
                (n_fets if nmos else p_fets)[anchor_node] += 1
                # NFet anchors its drain at the origin; for the PMOS we shift
                # left so the drain tip still lands on the node column
                base = (col + _H * index, 0.0) if nmos else (col - 1.5 - _H * index, 0.0)
                device = drawing.add(fet_cls().at(base).theta(90).label(label, loc="bot"))
                anchors = device.absanchors

                # drain: on the node dot for the first NMOS slot; raised lane otherwise
                if anchor_node == "0" or (nmos and index == 0):
                    dot(anchor_node)
                else:
                    _route(drawing, sym, anchors["drain"], col, _BODY_LANE)

                for role, node_name in _pin_nodes(element):
                    anchor_point = anchors.get(role)
                    if node_name == "0" and role != "gate":
                        drawing.add(sym.Ground().at(anchor_point))
                        continue
                    target_x = positions.get(node_name, (None,))[0]
                    if target_x is None or anchor_point is None:
                        continue
                    if role == "gate":
                        level = _GATE_LANE + _GATE_STEP * gate_lanes
                        gate_lanes += 1
                        _route(drawing, sym, anchor_point, target_x, level)
                    elif node_name != anchor_node:
                        level = _WIRE_LANE - _WIRE_STEP * wire_lanes
                        wire_lanes += 1
                        _route(drawing, sym, anchor_point, target_x, level)
                    dot(node_name)
                continue

            # ---------------- two-terminal elements ----------------
            symbol_cls = getattr(sym, _symbol_name(element) or "Resistor")
            n1, n2 = element.n1, element.n2
            p1, p2 = positions.get(n1), positions.get(n2)

            if n1 != "0" and n2 != "0" and p1 is not None and p2 is not None:
                device = drawing.add(symbol_cls().at(p1).to(p2).label(label, loc="top"))
                drawing.add(sym.Dot().at(device.absanchors["start"]))
                drawing.add(sym.Dot().at(device.absanchors["end"]))
                for node_name in (n1, n2):
                    if node_name not in labeled:
                        labeled.add(node_name)
                        drawing.add(
                            sym.Label().at(positions[node_name])
                            .label(node_name, loc="bot", fontsize=9, color="#555")
                        )
            elif n1 == "0" or n2 == "0":
                real_node = n1 if n1 != "0" else n2
                if real_node not in positions:
                    continue
                shift = _H * hangs[real_node]
                hangs[real_node] += 1
                start = (positions[real_node][0] + shift, 0.0)
                device = drawing.add(symbol_cls().at(start).down(_HANG).label(label, loc="bot"))
                if shift:
                    drawing.add(sym.Line().at(positions[real_node]).to(start))
                drawing.add(sym.Ground().at(device.absanchors["end"]))
                dot(real_node)
            # both pins on ground: nothing to draw

        drawing.save(path, dpi=110)
    return True
