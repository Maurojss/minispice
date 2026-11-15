---
name: minispice-new-device
description: Step-by-step recipe for adding a device model, element or analysis to MiniSPICE - parser entry, MNA stamp + Jacobian, detailed .op report, per-element current, transient/AC support and the mandatory closed-form validation test. Use when adding or extending devices (JFET, graded junction caps, CJS, flicker noise) or analyses.
---

# Adding a device or analysis to MiniSPICE

Project: `Projetos/minispice`. Read `README.md` ("How it works", "Limitations") first. Core files:
`src/minispice/parser.py` (syntax/params), `src/minispice/solver.py` (stamps, Newton, analyses),
`src/minispice/results.py` (result containers), `tests/`.

## The iron rule

**No feature without a closed-form reference test.** Pick a tiny circuit whose answer you can compute
analytically (or from the PDK/textbook equation), assert within a stated tolerance, and put the
reference value in the comment. Follow `tests/test_ac.py` / `tests/test_advanced.py` patterns.

## Recipe (new element + model)

1. **Parse** (`parser.py`):
   - Add the letter to the element dispatch; parse name, nodes, value/model, `TOL=`, `IC=`.
   - Add a model dataclass if needed (see `DiodeModel`, and how `NSUB` hangs off the MOSFET model).
   - Error messages: point at the offending line and the expectation.
2. **Stamp + Jacobian** (`solver.py`, inside `_assemble` / a `_stamp_*` method):
   - Linear part: stamp conductances into `A`.
   - Nonlinear part: linearize `I(v) ~= I(v0) + J(v0)(v - v0)`; `J` goes into `A`,
     `J v0 - I(v0)` into `b`.
   - Keep DC / transient / AC behavior consistent; in AC a device is its resistive Jacobian piece
     (or complex stamp) - decide deliberately and document.
3. **Charge storage** (if any): follow `_stamp_cap_branch` companion models.
   Capacitor: BE `G=C/dt, Ieq=-G v_prev`; trap `G=2C/dt, Ieq=-G v_prev - i_prev`. The transient
   companion model owns the stored charge.
4. **Per-element current**: register `I:<name>` so it shows in CLI/GUI/CSV (see `_diode_currents`,
   `_bjt_currents`). Reverse leakage must be exact where physics says so (diode: exactly `-IS`).
5. **Detailed `.op`**: extend `_device_details` - region, currents, small-signal quantities
   (`gm`, `gds`, `beta`, `Vth`, `Vov`, ...). These feed the GUI and the CLI table.
6. **Hidden internal terminals** for parasitic resistances: `_unique_node` + `<name>#<pin>` pattern
   (see `_build_bjt_parasitics`, `q1#b`); scale by area; never appear in the public node list.
7. **Tests** (one file or the thematic file):
   - Closed-form DC/transient/AC reference assertions.
   - A failure-path test (bad model param -> clear parse error).
   - If it changes a documented number: update the README *Validation* table in the same PR.
8. **Gates**:
   ```bash
   pytest && pytest --cov=minispice --cov-fail-under=75 && ruff check src tests scripts tools
   ```

## Recipe (new analysis)

- Add a method on `Simulation` (see `dc_sweep`, `sensitivity`), wire the dot-command in the parser and
  the CLI/GUI paths, return a typed result container from `results.py` with plot support.
- Validate against an analytic circuit (e.g. a new `.pz` idea: single-pole RC, compare to
  `1/(1+sRC)`).
- Update README (Features + example + Limitations) and the roadmap checkbox if it closes one.

## Things deliberately NOT in the project (do not add silently)

No BSIM charge model for the Level-1 MOSFET (transient capacitance handled externally via
`docs/capacitance.md`), no JFET or transmission lines yet, no flicker noise, no nested `.dc` sweeps,
no per-device MC distributions or Latin hypercube. Extending any of these = update *Limitations*.
