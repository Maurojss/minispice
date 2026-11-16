---
name: minispice-dev-workflow
description: Dev environment, quality gates and conventions for the MiniSPICE project - venv setup, pytest/coverage/ruff commands, code style, project layout, GUI run script and what must stay green before finishing any task.
---

# MiniSPICE dev workflow

Project: `Projetos/minispice` (Python 3.10+, NumPy core, optional matplotlib/schemdraw/scipy).

## Environment

```bash
cd Projetos/minispice
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pre-commit install            # optional: ruff hook on commit
```

`./gui.sh` prepares `.venv` if needed and launches the browser GUI (`minispice --gui` equivalent).

## Quality gates (run before reporting done)

```bash
pytest                                        # suite (currently 243 tests), addopts already -q
pytest --cov=minispice --cov-fail-under=75    # CI enforces 75% coverage (currently ~86%)
ruff check src tests scripts tools            # select = ["E9", "F", "I"], line-length 100
```

All three green = task done. A task that adds code without tests is not done.

## Layout

```text
src/minispice/   parser.py (~1.5k lines) solver.py (~2.3k) results.py schematic.py cli.py webui.py
cells/sky130/    transistor-level compound cells (inv nand2 nor2 aoi21 oai21 aoi22 oai22 mux21)
pdk/sky130/      vendored real sky130 BSIM4 files + extract.py + summary.md (calibration mapping)
examples/        runnable .cir netlists (one per feature, referenced from README)
scripts/         verify_gates.py (truth tables), aoi21_vtc.py (VTC plot)
tests/           pytest suite; closed-form reference pattern in test_ac.py / test_advanced.py
tools/           testbench helpers shared by tests and scripts
docs/            generated plots + worked notes (capacitance.md)
```

## Conventions

- **Didactic readability over cleverness** - every numerical step should be visible and documented.
- **Closed-form references in tests** with the analytic value cited in the assertion comment.
- **README is the spec**: features, validation table numbers, limitations and roadmap stay in sync
  with code in the same change.
- **Per-element currents as `I:<name>`** everywhere (CLI, GUI, CSV, Python).
- Netlist examples in `examples/` are runnable and referenced from README; new syntax ships with one.
- Commit style: keep commits scoped to one concern; run the gates before committing (pre-commit runs
  ruff anyway).

## CI

Runs the whole pytest suite on every push + coverage gate + ruff. Treat any CI red as blocking.
