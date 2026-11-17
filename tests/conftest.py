import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # so that `import tools.gates` works

from minispice import NetlistParser, Simulation  # noqa: E402


@pytest.fixture
def run():
    def _run(text: str):
        return Simulation(NetlistParser(text).parse()).run()

    return _run


@pytest.fixture
def eval_gate():
    """Run a compound cell for one input combination, return (V(y), logic).

    Uses fresh parsers per combination so every testbench is independent.
    """

    def _eval(gate, inputs):
        from tools.gates import build_testbench

        circuit = NetlistParser(build_testbench(gate, inputs)).parse()
        result = Simulation(circuit).operating_point()
        return result.nodes[gate.output]

    return _eval
