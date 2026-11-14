"""Command-line interface."""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

from . import __version__
from .parser import NetlistError, NetlistParser, SensitivitySpec
from .results import AcResult, DcSweep, NoiseResult, OperatingPoint, Transient
from .solver import ConvergenceError, Simulation


def _write_csv(result, path: Path) -> None:
    if isinstance(result, OperatingPoint):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["signal", "value"])
            for name, value in result.nodes.items():
                writer.writerow([f"V({name})", value])
            for name, value in result.currents.items():
                writer.writerow([f"I({name})", value])
    elif isinstance(result, DcSweep):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow([result.source, *[f"V({n})" for n in result.nodes], *[f"I({n})" for n in result.currents]])
            for i, value in enumerate(result.sweep):
                writer.writerow([value, *[series[i] for series in result.nodes.values()], *[series[i] for series in result.currents.values()]])
    elif isinstance(result, Transient):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(["time", *[f"V({n})" for n in result.nodes], *[f"I({n})" for n in result.currents]])
            for i, value in enumerate(result.time):
                writer.writerow([value, *[series[i] for series in result.nodes.values()], *[series[i] for series in result.currents.values()]])
    elif isinstance(result, AcResult):
        pairs = [
            *((f"V({name})", series) for name, series in result.nodes.items()),
            *((f"I({name})", series) for name, series in result.currents.items()),
        ]
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(
                [
                    "freq",
                    *(f"|{label}| dB" for label, _ in pairs),
                    *(f"arg({label}) deg" for label, _ in pairs),
                ]
            )
            for i, value in enumerate(result.freqs):
                mags = [
                    20.0 * math.log10(max(abs(series[i]), 1e-30))
                    for _, series in pairs
                ]
                phases = [math.degrees(math.atan2(series[i].imag, series[i].real))
                          for _, series in pairs]
                writer.writerow([value, *mags, *phases])
    elif isinstance(result, NoiseResult):
        has_input = result.input_psd is not None
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            header = ["freq", "output_psd", "output_density"]
            if has_input:
                header += ["input_psd", "input_density"]
            writer.writerow(header)
            for index, (value, psd) in enumerate(
                zip(result.freqs, result.output_psd)
            ):
                row = [value, psd, math.sqrt(max(psd, 0.0))]
                if has_input and result.input_psd is not None:
                    in_psd = result.input_psd[index]
                    row += [in_psd, math.sqrt(max(in_psd, 0.0))]
                writer.writerow(row)


def _write_summary_csv(summary: dict[str, dict[str, float]], path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["signal", "mean", "std", "min", "max"])
        for name, stats in summary.items():
            writer.writerow(
                [name, stats["mean"], stats["std"], stats["min"], stats["max"]]
            )


def _print_sensitivity(reports: list[dict]) -> None:
    for report in reports:
        print(f"Sensitivity: {report['output']}")
        print(f"  nominal = {report['nominal']:.6g}")
        print(f"  {'device':<12}{'value':>12}{'dY/dP':>14}{'normalized':>14}")
        for entry in report["entries"]:
            print(
                f"  {entry['name']:<12}{entry['value']:>12.6g}"
                f"{entry['derivative']:>14.4g}{entry['sensitivity']:>14.4g}"
            )
        worst = report.get("worst_case") or {}
        low, high = worst.get("min"), worst.get("max")
        if low and high:
            def label(point: dict) -> str:
                return ", ".join(
                    f"{name}={value:g}" for name, value in point["case"].items()
                ) or "nominal"

            print(
                f"  worst case: min {low['value']:.6g} ({label(low)})"
                f" / max {high['value']:.6g} ({label(high)})"
            )
        print()


def _plot_path(base: Path, index: int, total: int) -> Path:
    if total == 1:
        return base
    return base.with_name(f"{base.stem}_{index}{base.suffix or '.png'}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="minispice",
        description="Small SPICE-like circuit simulator (MNA + Newton-Raphson).",
    )
    parser.add_argument("netlist", nargs="?", help="path to the .cir netlist")
    parser.add_argument(
        "--version", action="version", version=f"minispice {__version__}"
    )
    parser.add_argument("--gui", action="store_true",
                        help="open the interactive browser GUI (netlist editor + plots)")
    parser.add_argument("--no-browser", action="store_true",
                        help="with --gui: do not open the browser automatically")
    parser.add_argument("--port", type=int, default=None, metavar="N",
                        help="GUI server port (default 8080)")
    parser.add_argument("--plot", metavar="PATH", help="save plots to PATH (PNG)")
    parser.add_argument("--show", action="store_true", help="show plots interactively")
    parser.add_argument("--csv", metavar="PATH",
                        help="export results to CSV (Monte Carlo writes one summary CSV)")
    parser.add_argument(
        "--sensitivity",
        metavar="SIGNAL",
        help="rank DC operating-point sensitivities and report worst-case corners",
    )
    parser.add_argument(
        "--worst-case",
        metavar="SIGNAL",
        help="alias for --sensitivity (the report always includes worst-case corners)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="only errors")
    args = parser.parse_args(argv)

    if args.gui:
        from . import webui

        gui_args: list[str] = []
        if args.port is not None:
            gui_args += ["--port", str(args.port)]
        if args.no_browser:
            gui_args.append("--no-browser")
        if args.netlist:
            gui_args.append(args.netlist)
        return webui.main(gui_args)

    try:
        if not args.netlist:
            parser.error("the following arguments are required: netlist")
        circuit = NetlistParser.from_file(args.netlist).parse()
        signal = args.sensitivity or args.worst_case
        if signal:
            circuit.sensitivity.append(SensitivitySpec(signal))
        simulation = Simulation(circuit)
        results = simulation.run()
    except (NetlistError, ConvergenceError, ValueError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    if not args.quiet:
        if circuit.title:
            print(circuit.title)
            print("=" * len(circuit.title))
        for warning in circuit.warnings:
            print(f"warning: {warning}", file=sys.stderr)

    total = len(results)
    summary = getattr(simulation, "mc_summary", {})
    if summary and not args.quiet:
        print(f"Monte Carlo summary ({circuit.mc_runs} runs):")
        print(f"  {'signal':<24}{'mean':>12}{'std':>12}{'min':>12}{'max':>12}")
        for name, stats in summary.items():
            print(
                f"  {name:<24}{stats['mean']:>12.6g}{stats['std']:>12.6g}"
                f"{stats['min']:>12.6g}{stats['max']:>12.6g}"
            )
        print()

    yields = getattr(simulation, "yield_summary", {})
    if yields and not args.quiet:
        print("Monte Carlo yield:")
        for name, stats in yields.items():
            if stats.get("yield") is None:
                continue
            bounds = []
            if stats.get("limit_min") is not None:
                bounds.append(f"min {stats['limit_min']:g}")
            if stats.get("limit_max") is not None:
                bounds.append(f"max {stats['limit_max']:g}")
            print(
                f"  {name:<24}{stats['yield'] * 100:>6.1f}%"
                f" ({stats['pass']}/{stats['total']})"
                f"  {' '.join(bounds)}"
            )
        print()

    sensitivity = getattr(simulation, "sensitivity_report", [])
    if sensitivity and not args.quiet:
        _print_sensitivity(sensitivity)

    if args.csv and summary:
        path = Path(args.csv).with_suffix(".csv")
        path.parent.mkdir(parents=True, exist_ok=True)
        _write_summary_csv(summary, path)
        if not args.quiet:
            print(f"csv:  {path}")
            print()

    for index, result in enumerate(results, start=1):
        if not args.quiet and not summary:
            print(result)
            print()

        if args.plot and not isinstance(result, OperatingPoint):
            path = _plot_path(Path(args.plot), index, total)
            path.parent.mkdir(parents=True, exist_ok=True)
            result.plot(path=str(path), show=False)
            if not args.quiet:
                print(f"plot: {path}")

        if args.csv and not summary:
            path = _plot_path(Path(args.csv), index, total).with_suffix(".csv")
            path.parent.mkdir(parents=True, exist_ok=True)
            _write_csv(result, path)
            if not args.quiet:
                print(f"csv:  {path}")

    return 0
