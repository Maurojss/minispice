"""Simulation result containers and plotting helpers."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _format_measurements(measurements: dict[str, float]) -> list[str]:
    if not measurements:
        return []
    lines = ["  measurements:"]
    for name, value in measurements.items():
        shown = "failed" if not np.isfinite(value) else f"{value:.6g}"
        lines.append(f"    {name} = {shown}")
    return lines


def _format_case(case: dict[str, float]) -> str:
    return "[" + ", ".join(f"{name}={value:g}" for name, value in case.items()) + "]"


@dataclass
class OperatingPoint:
    nodes: dict[str, float]
    currents: dict[str, float]
    iterations: int = 0
    vector: np.ndarray | None = field(default=None, repr=False)
    details: dict[str, dict[str, float | str]] = field(default_factory=dict)
    measurements: dict[str, float] = field(default_factory=dict)
    case: dict[str, float] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = ["Operating point:"]
        if self.case:
            lines[0] += " " + _format_case(self.case)
        for name, value in self.nodes.items():
            lines.append(f"  V({name}) = {value:.6g} V")
        for name, value in self.currents.items():
            lines.append(f"  I({name}) = {value:.6g} A")
        if self.details:
            lines.append("  devices:")
            for name, info in self.details.items():
                region = info.get("region", "")
                if info.get("type") in ("nmos", "pmos"):
                    lines.append(
                        f"    {name} ({info.get('type')}, {region}): "
                        f"Id={float(info['id']):.6g} A  "
                        f"Vgs={float(info['vgs']):.4g} V  "
                        f"Vds={float(info['vds']):.4g} V  "
                        f"Vth={float(info['vth']):.4g} V  "
                        f"gm={float(info['gm']):.4g} S  "
                        f"gds={float(info['gds']):.4g} S"
                    )
                elif info.get("type") == "diode":
                    lines.append(
                        f"    {name} (diode, {region}): "
                        f"Vd={float(info['vd']):.6g} V  "
                        f"Id={float(info['id']):.6g} A  "
                        f"gd={float(info['gd']):.4g} S"
                    )
                else:
                    lines.append(
                        f"    {name} ({info.get('type')}, {region}): "
                        f"Ic={float(info['ic']):.6g} A  "
                        f"Ib={float(info['ib']):.6g} A  "
                        f"Vbe={float(info['vbe']):.4g} V  "
                        f"Vce={float(info['vce']):.4g} V  "
                        f"beta={float(info['beta']):.4g}"
                    )
        lines.extend(_format_measurements(self.measurements))
        lines.append(f"  iterations: {self.iterations}")
        return "\n".join(lines)


@dataclass
class DcSweep:
    source: str
    sweep: np.ndarray
    nodes: dict[str, np.ndarray]
    currents: dict[str, np.ndarray]
    iterations: list[int] = field(default_factory=list)
    measurements: dict[str, float] = field(default_factory=dict)
    case: dict[str, float] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [f"DC sweep of {self.source}:"]
        if self.case:
            lines[0] += " " + _format_case(self.case)
        header = f"{self.source:>12}"
        for name in self.nodes:
            header += f"  V({name})"
        lines.append(header)
        for i, value in enumerate(self.sweep):
            row = f"{value:12.6g}"
            for series in self.nodes.values():
                row += f"  {series[i]:12.6g}"
            lines.append(row)
        lines.extend(_format_measurements(self.measurements))
        return "\n".join(lines)

    def plot(self, nodes: list[str] | None = None, path: str | None = None, show: bool = False):
        return _plot_series(self.sweep, self.nodes, self.source, "Voltage (V)", nodes, path, show)


@dataclass
class Transient:
    time: np.ndarray
    nodes: dict[str, np.ndarray]
    currents: dict[str, np.ndarray]
    iterations: list[int] = field(default_factory=list)
    measurements: dict[str, float] = field(default_factory=dict)
    case: dict[str, float] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [
            f"Transient: {len(self.time)} points, t = {self.time[0]:.6g} .. {self.time[-1]:.6g} s"
        ]
        if self.case:
            lines[0] += " " + _format_case(self.case)
        if self.nodes:
            header = f"{'time':>12}"
            for name in self.nodes:
                header += f"  V({name})"
            lines.append(header)
            for idx in (0, len(self.time) - 1):
                row = f"{self.time[idx]:12.6g}"
                for series in self.nodes.values():
                    row += f"  {series[idx]:12.6g}"
                lines.append(row)
        lines.extend(_format_measurements(self.measurements))
        return "\n".join(lines)

    def plot(self, nodes: list[str] | None = None, path: str | None = None, show: bool = False):
        return _plot_series(self.time, self.nodes, "Time (s)", "Voltage (V)", nodes, path, show)

    def _series(self, signal: str) -> np.ndarray:
        text = signal.strip()
        lowered = text.lower()
        if lowered.startswith("v(") and lowered.endswith(")"):
            key = text[2:-1]
            return np.asarray(self.nodes[key], dtype=float)
        if lowered.startswith("i(") and lowered.endswith(")"):
            key = text[2:-1]
            return np.asarray(self.currents[key], dtype=float)
        if text in self.nodes:
            return np.asarray(self.nodes[text], dtype=float)
        return np.asarray(self.currents[text], dtype=float)

    def fft(self, signal: str = "out") -> tuple[np.ndarray, np.ndarray]:
        """Single-sided amplitude spectrum of a transient signal."""
        y = self._series(signal)
        n = len(y)
        if n < 2:
            raise ValueError("FFT needs at least two samples")
        dt = float(np.mean(np.diff(self.time)))
        if dt <= 0:
            raise ValueError("time vector must be increasing")
        window = y - float(np.mean(y))
        spectrum = np.abs(np.fft.rfft(window)) * (2.0 / n)
        freqs = np.fft.rfftfreq(n, d=dt)
        return freqs, spectrum

    def thd(self, frequency: float, signal: str = "out", harmonics: int = 5) -> float:
        """Total harmonic distortion from harmonics 2..N of ``frequency``."""
        freqs, spectrum = self.fft(signal)
        dt = float(np.mean(np.diff(self.time)))
        resolution = 1.0 / (len(self.time) * dt)

        def amplitude(order: int) -> float:
            target = order * frequency
            mask = np.abs(freqs - target) <= resolution
            return float(np.max(spectrum[mask])) if mask.any() else 0.0

        fundamental = amplitude(1)
        if fundamental <= 0.0:
            return float("nan")
        return float(np.sqrt(sum(amplitude(k) ** 2 for k in range(2, harmonics + 1))) / fundamental)


@dataclass
class AcResult:
    """Small-signal frequency sweep: complex phasors at each frequency."""

    freqs: np.ndarray
    nodes: dict[str, np.ndarray]      # complex voltage phasors
    currents: dict[str, np.ndarray]   # complex branch-current phasors
    measurements: dict[str, float] = field(default_factory=dict)
    case: dict[str, float] = field(default_factory=dict)

    def db(self, series: np.ndarray) -> np.ndarray:
        """Magnitude in dB (0 dB reference used instead of dBV)."""
        magnitude = np.abs(series)
        return 20.0 * np.log10(np.where(magnitude > 0, magnitude, np.nan))

    def phase(self, series: np.ndarray) -> np.ndarray:
        """Phase in degrees (``arg`` of the phasor)."""
        return np.rad2deg(np.angle(series))

    def __str__(self) -> str:
        lines = [f"AC analysis: {len(self.freqs)} points, "
                 f"f = {self.freqs[0]:.6g} .. {self.freqs[-1]:.6g} Hz"]
        if self.case:
            lines[0] += " " + _format_case(self.case)
        if self.nodes:
            header = "  freq(Hz)".rjust(12)
            for name in self.nodes:
                header += f"  |V({name})|dB"
            lines.append(header)
            for idx in (0, len(self.freqs) - 1):
                row = f"{self.freqs[idx]:12.6g}"
                for series in self.nodes.values():
                    mag = abs(series[idx])
                    row += f"  {(20 * np.log10(mag) if mag > 0 else float('-inf')):12.4g}"
                lines.append(row)
        lines.extend(_format_measurements(self.measurements))
        return "\n".join(lines)

    def plot(self, nodes: list[str] | None = None, path: str | None = None, show: bool = False):
        import matplotlib

        if path or not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        names = nodes or list(self.nodes)
        fig, ax = plt.subplots(figsize=(8, 4.5))
        for name in names:
            if name not in self.nodes:
                raise KeyError(f"unknown signal: {name}")
            ax.plot(self.freqs, self.db(self.nodes[name]), label=name)
        ax.set_xscale("log")
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Magnitude (dB)")
        ax.grid(True, which="both", alpha=0.3)
        if names:
            ax.legend()
        fig.tight_layout()
        if path:
            fig.savefig(path, dpi=140)
            plt.close(fig)
        if show:
            plt.show()
        return fig


@dataclass
class NoiseResult:
    """Output noise power spectral density (V^2/Hz)."""

    freqs: np.ndarray
    output_psd: np.ndarray
    output_name: str = "out"
    input_source: str = ""
    input_unit: str = "V"
    total_rms: float = 0.0
    input_psd: np.ndarray | None = None
    input_rms: float = 0.0
    measurements: dict[str, float] = field(default_factory=dict)
    case: dict[str, float] = field(default_factory=dict)

    def __str__(self) -> str:
        lines = [
            f"Noise analysis at V({self.output_name}): {len(self.freqs)} points, "
            f"f = {self.freqs[0]:.6g} .. {self.freqs[-1]:.6g} Hz"
        ]
        if self.case:
            lines[0] += " " + _format_case(self.case)
        lines.append(f"  integrated output noise: {self.total_rms:.6g} V rms")
        if self.input_psd is not None:
            lines.append(
                f"  integrated input noise:  {self.input_rms:.6g} "
                f"{self.input_unit} rms"
            )
        for name, value in self.measurements.items():
            lines.append(f"  {name} = {value:.6g}")
        return "\n".join(lines)

    def plot(self, path: str | None = None, show: bool = False):
        import matplotlib

        if path or not show:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4.5))
        psd = np.where(self.output_psd > 0, self.output_psd, np.nan)
        ax.loglog(self.freqs, np.sqrt(psd), label="output noise")
        if self.input_psd is not None:
            in_psd = np.where(self.input_psd > 0, self.input_psd, np.nan)
            ax.loglog(
                self.freqs,
                np.sqrt(in_psd),
                label=f"input-referred noise ({self.input_unit}/sqrt(Hz))",
            )
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Noise density (V/sqrt(Hz))")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        if path:
            fig.savefig(path, dpi=140)
            plt.close(fig)
        if show:
            plt.show()
        return fig


def _plot_series(
    x: np.ndarray,
    series: dict[str, np.ndarray],
    xlabel: str,
    ylabel: str,
    nodes: list[str] | None,
    path: str | None,
    show: bool,
):
    import matplotlib

    if path:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    names = nodes or list(series)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for name in names:
        if name not in series:
            raise KeyError(f"unknown signal: {name}")
        ax.plot(x, series[name], label=name)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    if names:
        ax.legend()
    fig.tight_layout()
    if path:
        fig.savefig(path, dpi=140)
        plt.close(fig)
    if show:
        plt.show()
    return fig
