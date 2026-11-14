# Contributing to MiniSPICE

Thanks for your interest! MiniSPICE is first of all a **learning project**:
readable code, visible numerics and tests that pin every result to a
closed-form value. Those three things make a good contribution.

## Development setup

```bash
git clone https://github.com/Maurojss/minispice.git
cd minispice
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

The `dev` extra pulls pytest plus matplotlib, schemdraw and esprima (used to
validate the GUI page JavaScript from inside the test suite).

## Running the checks

```bash
pytest                       # full suite (~129 tests, should stay green)
pytest tests/test_ac.py -q   # just one area
./gui.sh                     # open the browser GUI; prepares .venv if needed
./gui.sh --port 8081         # use an alternate port
# Or, after installation: minispice --gui
```

If the GUI behaves differently from what you see, a screenshot in the issue
helps a lot.

## Project conventions

- **Style:** plain, documented Python. Docstrings on every public function,
  comments explaining *why* a numerical trick is used (companion models, gmin,
  damping...). Line length is soft ~100 chars; imports at the top of the file.
- **Modules:**
  - `parser.py` — text -> `Circuit` (no numerics)
  - `solver.py` — all numerics (MNA, Newton, transient, AC)
  - `results.py` — result containers and their `__str__`/`plot()`
  - `webui.py` / `schematic.py` — the browser GUI; keep it dependency-light
    (stdlib server + matplotlib + schemdraw)
- **Netlists in tests are inline strings**, so each test is self-contained
  like `LOWPASS = """..."""` in `tests/test_ac.py`.

## Adding a numerical feature

Every solver change should come with an analytic reference:

1. Write the closed-form answer first (e.g. `|H(fc)| = -3.0103 dB` in the AC
   tests, `1 - e^-1` in the RC step test).
2. Assert the simulator against it with a defensible tolerance.
3. Mention the validation basis in the PR description — this is what makes
   the result trustworthy rather than just plausible.

For new device models, follow the pattern of `_stamp_diode`/`_stamp_mosfet`
in `src/minispice/solver.py`: a `_xxx_params` helper that returns the
linearized parameters, reused by both DC (`_stamp_nonlinear`) and AC
(small-signal stamps).

## Commits and PRs

- Small, focused PRs. One behavioral change per PR.
- Include the commands you used to verify (`pytest`, `minispice ...`).
- Screenshots or PNGs are welcome for plot/GUI changes.

## Reporting bugs

Open an issue with:

- the **netlist** (several lines, runnable, please paste it);
- the command you ran and the full output;
- the Python version: `python --version`;
- what you expected vs. what happened.

License of contributions: MIT (same as the project).
