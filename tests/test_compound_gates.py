import pytest

from tools.gates import GATES, HIGH, LOW, VDD, build_testbench, combos

ALL_GATES = list(GATES.values())
# One (gate, inputs) case per truth-table row: 2+4+4+8+8+16+16+8 = 66 checks.
ALL_CASES = [(g, c) for g in ALL_GATES for c in combos(g)]
CASE_IDS = [f"{g.name}-{''.join(map(str, c))}" for g, c in ALL_CASES]


def test_registry_matches_netlist_ports():
    """Declared inputs and output `y` must exist in the flattened circuit."""
    from minispice import NetlistParser

    for gate in ALL_GATES:
        circuit = NetlistParser(
            build_testbench(gate, (0,) * len(gate.inputs))
        ).parse()
        assert "y" in circuit.node_names
        for port in gate.inputs:
            assert port in circuit.node_names, (gate.name, port)


def test_margins_are_valid():
    assert 0.0 <= LOW < 0.5 * VDD < HIGH <= VDD


def level_of(y: float) -> int | None:
    if y <= LOW:
        return 0
    if y >= HIGH:
        return 1
    return None


@pytest.mark.parametrize(("gate", "inputs"), ALL_CASES, ids=CASE_IDS)
def test_truth_table(eval_gate, gate, inputs):
    y = eval_gate(gate, inputs)
    got = level_of(y)
    assert got is not None, f"{gate.name}{inputs}: ambiguous y={y:.4g} V"
    assert got == gate.expect(*inputs), (
        f"{gate.name}{inputs}: y={y:.4g} V -> {got}, "
        f"expected {gate.expect(*inputs)}"
    )
