"""MiniSPICE — a small, readable SPICE-like circuit simulator."""

from .parser import NetlistError, NetlistParser
from .results import AcResult, DcSweep, NoiseResult, OperatingPoint, Transient
from .solver import ConvergenceError, Simulation

__all__ = [
    "NetlistParser",
    "NetlistError",
    "Simulation",
    "ConvergenceError",
    "OperatingPoint",
    "DcSweep",
    "Transient",
    "AcResult",
    "NoiseResult",
]

__version__ = "0.1.0"
