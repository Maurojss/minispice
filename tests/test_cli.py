import csv
import runpy
import subprocess
import sys

import pytest

from minispice.cli import main


def test_cli_csv(tmp_path):
    netlist = tmp_path / "rc.cir"
    netlist.write_text(
        "V1 in 0 PULSE(0 1 0 1n 1n 100m 200m)\n"
        "R1 in out 1k\nC1 out 0 1u\n.tran 10u 1m\n.end\n"
    )
    out = tmp_path / "result.csv"
    assert main([str(netlist), "--csv", str(out), "--quiet"]) == 0
    assert out.exists()
    with out.open() as fh:
        rows = list(csv.reader(fh))
    assert rows[0][0] == "time"
    assert len(rows) > 10


def test_cli_plot(tmp_path):
    netlist = tmp_path / "divider.cir"
    netlist.write_text("V1 1 0 10\nR1 1 2 1k\nR2 2 0 1k\n.op\n.end\n")
    out = tmp_path / "divider.png"
    assert main([str(netlist), "--plot", str(out), "--quiet"]) == 0
    assert not out.exists()  # .op has no plot


def test_cli_plot_transient(tmp_path):
    netlist = tmp_path / "rc.cir"
    netlist.write_text(
        "V1 in 0 PULSE(0 1 0 1n 1n 10m 20m)\n"
        "R1 in out 1k\nC1 out 0 1u\n.tran 10u 1m\n.end\n"
    )
    out = tmp_path / "rc.png"
    assert main([str(netlist), "--plot", str(out), "--quiet"]) == 0
    assert out.exists() and out.stat().st_size > 0


def test_cli_noise_csv_and_plot(tmp_path):
    netlist = tmp_path / "noise.cir"
    netlist.write_text(
        "V1 in 0 dc 0 ac 1\n"
        "R1 in out 1k\nC1 out 0 1u\n"
        ".noise V(out) V1 dec 5 1 1e6\n.end\n"
    )
    csv_out = tmp_path / "noise.csv"
    assert main([str(netlist), "--csv", str(csv_out), "--quiet"]) == 0
    with csv_out.open() as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == [
        "freq",
        "output_psd",
        "output_density",
        "input_psd",
        "input_density",
    ]
    assert len(rows) == 32

    png_out = tmp_path / "noise.png"
    assert main([str(netlist), "--plot", str(png_out), "--quiet"]) == 0
    assert png_out.exists() and png_out.stat().st_size > 0


def test_cli_monte_carlo_summary(tmp_path, capsys):
    netlist = tmp_path / "mc.cir"
    netlist.write_text(
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".dc V1 1 1 1\n"
        ".mc 5 SEED=2\n.end\n"
    )
    assert main([str(netlist)]) == 0
    output = capsys.readouterr().out
    assert "Monte Carlo summary (5 runs)" in output
    assert "V(out)" in output

    csv_out = tmp_path / "mc_summary.csv"
    assert main([str(netlist), "--csv", str(csv_out), "--quiet"]) == 0
    with csv_out.open() as fh:
        rows = list(csv.reader(fh))
    assert rows[0] == ["signal", "mean", "std", "min", "max"]
    assert any(row[0] == "V(out)" for row in rows[1:])
    assert not (tmp_path / "mc_summary_1.csv").exists()


def test_cli_reports_monte_carlo_yield(tmp_path, capsys):
    netlist = tmp_path / "yield.cir"
    netlist.write_text(
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".dc V1 0.5 1.5 0.5\n"
        ".meas dc vout MAX V(out) MIN=0.4 MAX=0.9\n"
        ".mc 5 SEED=1\n.end\n"
    )
    assert main([str(netlist)]) == 0
    output = capsys.readouterr().out
    assert "Monte Carlo yield:" in output
    assert "vout" in output and "100.0%" in output


def test_cli_sensitivity_flag(tmp_path, capsys):
    netlist = tmp_path / "sens.cir"
    netlist.write_text(
        "Divider\n"
        "V1 in 0 1\n"
        "R1 in out 1k TOL=0.1\n"
        "R2 out 0 1k TOL=0.1\n"
        ".op\n.end\n"
    )
    assert main([str(netlist), "--sensitivity", "V(out)"]) == 0
    output = capsys.readouterr().out
    assert "Sensitivity: V(out)" in output
    assert "R1" in output and "R2" in output
    assert "worst case:" in output
    assert "0.45" in output and "0.55" in output


def test_python_module_entry_point_help_and_netlist(tmp_path, monkeypatch, capsys):
    """``__main__.py`` run as a script: --help and a real netlist both exit 0.

    runpy executes ``src/minispice/__main__.py`` with ``__name__ == "__main__"``
    in this process, so the module-level import, the guard and
    ``raise SystemExit(main())`` are all exercised (a subprocess would leave
    them invisible to coverage). sys.argv is patched exactly as the CLI sees it.
    """
    netlist = tmp_path / "divider.cir"
    netlist.write_text("Divider\nV1 1 0 10\nR1 1 2 1k\nR2 2 0 1k\n.op\n.end\n")

    monkeypatch.setattr(sys, "argv", ["minispice", "--help"])
    with pytest.raises(SystemExit) as helper:
        runpy.run_module("minispice", run_name="__main__")
    assert helper.value.code == 0
    out = capsys.readouterr().out
    assert "usage: minispice" in out
    assert ".cir netlist" in out

    monkeypatch.setattr(sys, "argv", ["minispice", str(netlist)])
    with pytest.raises(SystemExit) as run:
        runpy.run_module("minispice", run_name="__main__")
    assert run.value.code == 0
    # Closed form: 10 V split by two equal resistors -> V(2) = 5 V.
    assert "V(2) = 5 V" in capsys.readouterr().out


def test_python_dash_m_subprocess_help(tmp_path):
    """`python -m minispice --help` works through the interpreter's -m machinery.

    The subprocess runs the same interpreter as the suite from a temporary cwd,
    proving the installed package (and not the repository layout) provides the
    entry point. --help exits through argparse's SystemExit(0), so both the
    return code and the usage banner are checked.
    """
    done = subprocess.run(
        [sys.executable, "-m", "minispice", "--help"],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=120,
    )
    assert done.returncode == 0
    assert "usage: minispice" in done.stdout
    assert ".cir netlist" in done.stdout
