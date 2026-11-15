---
name: minispice-netlist-ref
description: Compact MiniSPICE netlist syntax reference - elements, sources, models, dot commands (.op/.dc/.tran/.ac/.noise/.four/.param/.meas/.step/.mc/.sens), expressions and the Python API. Use whenever writing, reviewing or debugging a .cir file or calling the simulator from Python.
---

# MiniSPICE netlist reference

Project: `Projetos/minispice` (educational SPICE-like simulator, Python + NumPy).
Full docs: `README.md` in the project root. First line of a `.cir` file is always a comment.

## Elements

```spice
R1  n1  n2  1k  TOL=0.05        ; resistor, optional MC tolerance
C1  out 0   1u  IC=0.5          ; capacitor, initial condition
L1  a   b   1m                  ; inductor (branch current reported)
D1  out 0   dmod                ; diode -> .model
M1  d  g  s  b  NMOS W=1u L=1u  ; Level-1 MOSFET (d g s b order!)
Q1  c  b  0  qmod               ; NPN/PNP BJT (c b e order)
E1  out 0  c 0  1               ; VCVS  (V = gain * V(c)-V(0))
G1  out 0  in 0  1m              ; VCCS
F1  out 0  V1 2                 ; CCCS (control = voltage-source name)
H1  out 0  V1 1k                ; CCVS
B1  out 0  V={2*V(in)+0.5*V(in)*V(in)}   ; behavioral (V(node), V(a,b), I(el), time, temper)
X1  in  out  rclow              ; subcircuit instance
```

Value suffixes: `k`, `meg` (1e6, NOT `m`), `mil`, `u`, `n`, `p`, `f`, `G`, `T`.

## Sources

```spice
V1 in 0 dc 5              ; DC
V1 in 0 SIN(0 1 1k)       ; sin(offset amp freq)
V1 in 0 PULSE(0 1 0 10u 10u 2m 4m)  ; PULSE(v1 v2 delay rise fall width period)
V1 in 0 PWL(0 0 1m 1 2m 1)          ; piecewise linear
V1 in 0 dc 0 ac 1 0       ; + small-signal phasor (needed for .ac/.noise)
```

## Dot commands

```spice
.op                              ; detailed bias point (region, gm, beta, ...)
.dc VIN 0 1.8 0.01               ; DC sweep (no nesting; use .step for corners)
.tran 10u 5m                     ; step stop   (BE or trapezoidal)
.tran 10u 5m uic                 ; start from .ic instead of DC op
.ac dec 20 10 100K               ; dec|oct|lin points fstart fstop
.noise V(out) V1 dec 20 1 1e9    ; thermal + shot noise, input- and output-referred
.four 1k V(out)                  ; harmonics 1-5 + THD after a .tran
.temp 100                        ; analysis temperature
.ic V(out)=0.5                   ; initial conditions
.param R=1k C={2*R}
.func sq(x) = {x*x}
.include models.spice
.subckt rclow in out  ...  .ends
```

## Statistical / spec analyses

```spice
.step param R LIST 1k 2k 5k      ; .step temp 0 27 100 ; .step V1 0 5 1 (cartesian product)
.mc 200 SEED=7 DIST=GAUSS SIGMA=0.05 RHO=0.5   ; .options mc_tol=0.02 (uniform default 5%)
.meas tran trise TRIG v(out) VAL=0.1 RISE=1 TARG v(out) VAL=0.9 RISE=1
.meas ac gain GAIN V(out) AT=1k
.meas dc vout MAX V(out) MIN=0.45 MAX=0.55     ; spec -> yield summary under .mc
.sens V(out) TOP=8 DEVICES=6     ; DC sensitivities + worst-case corners (TOL= per device)
```

`.meas` kinds: `AVG MAX MIN RMS RISE FALL WIDTH CROSS GAIN PHASE BW LOOPGAIN`
(+`FROM TO VAL LOW HIGH AT DROP`).

## Expressions `{...}`

`+ - * / **`, parentheses, `sqrt exp log sin ...`, suffixes (`{2*1k}`), comparisons
(`< <= > >= == !=`), `and`/`or`/`not`, `if(cond, a, b)` (lazy), `gauss(mean, sigma)`
(deterministic seed; `.param`/`.func` only, never in runtime `B` expressions).

## Models

```spice
.model NMOS NMOS(VTO=0.5 KP=120u LAMBDA=0.02 GAMMA=0.4 PHI=0.7 NSUB=0)
.model PMOS PMOS(VTO=-0.5 KP=60u LAMBDA=0.02 GAMMA=0.4 PHI=0.7)
.model dmod D(IS=1e-14 N=1 RS=0 XTI=3 EG=1.11 TNOM=27)
.model qmod NPN(IS=1e-15 BF=100 VAF=100 IKF=0.15m ISE=1e-14 NE=1.5 RB=1k RC=20 RE=5 CJE=20p CJC=5p TF=2n TR=10n)
```

## CLI and Python

```bash
minispice examples/rc_lowpass.cir --plot rc.png     # --csv PATH --quiet --show
minispice net.cir --sensitivity V(out)              # or --worst-case V(out)
minispice --gui [--port 8081] [--no-browser]
```

```python
from minispice import NetlistParser, Simulation

circuit = NetlistParser.from_file("examples/rc_lowpass.cir").parse()
sim = Simulation(circuit)
r = sim.transient(10e-6, 5e-3)          # or sim.run() for all declared analyses
print(r.nodes["out"][-1], r.measurements)
r.plot(["in", "out"], path="rc.png")
freqs, spectrum = r.fft("v(out)"); thd = r.thd(1e3, "v(out)")
```

Result objects carry `case` (`.step` values or `{"mc": n}`) and `measurements`; MC runs populate
`sim.mc_summary` and `sim.yield_summary` / `yield_by_case`; `.sens` populates `sim.sensitivity_report`.

## Reporting currents without ammeters

Every element current is available as `I:<name>` in CLI, GUI, CSV and Python
(`I:V1` branch, `I:D1` diode - exact `-IS` in reverse, `I:Q1` BJT).
