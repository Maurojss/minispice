# OpenCode agents & skills (MiniSPICE)

Configuração local do projeto para o OpenCode (`opencode` CLI). Carrega automaticamente
quando uma sessão roda dentro de `Projetos/minispice`.

## Agents (`.opencode/agents/`)

Lançar com o tool `subagent` ou pedir: "use o X para ...".

| Agent | Área |
|---|---|
| `spice-solver-dev` | `solver.py`: MNA, Newton-Raphson, gmin/source stepping, companion models BE/trap, AC/noise/sensibilidade |
| `device-modeler` | Física dos dispositivos: diodo (Shockley/XTI), MOSFET Level-1 (NSUB), BJT (Ebers-Moll + GP subset), PDK sky130 |
| `netlist-reviewer` | `parser.py`: sintaxe SPICE, `.subckt`, `.param`/`.func`, expressões `{...}`, dot-commands, revisão de `.cir` |
| `testbench-engineer` | Suíte pytest, testes com valor analítico de referência, cobertura >= 75%, ruff |
| `cell-designer` | Células digitais transistor-level em `cells/sky130/` (registro em `tools/gates.py`, truth tables, README) |
| `cell-verification` | Verificação das células: truth tables completas, corners temp/VDD, Monte Carlo yield (`.meas MIN/MAX`), noise margins da VTC |
| `cell-char` | Caracterização liberty-like: VTC, CGG, delay `0.69·Req·CL` (método Weste & Harris, `docs/capacitance.md`), energia, corners |

## Skills (`.opencode/skills/`)

Carregar com o tool `skill` (ou mencionar `@minispice-netlist-ref` no prompt):

| Skill | Conteúdo |
|---|---|
| `minispice-netlist-ref` | Referência compacta de sintaxe de netlist + CLI + API Python |
| `minispice-new-device` | Receita para adicionar dispositivo/modelo/análise com teste de forma fechada |
| `minispice-dev-workflow` | Ambiente, gates de qualidade (pytest/cov/ruff), layout e convenções |

## Fluxo sugerido

### Núcleo do simulador
1. Ideia → `spice-solver-dev` (numérica) ou `device-modeler` (física) spec + implementação.
2. `testbench-engineer` escreve/roda os testes com referência analítica.
3. `netlist-reviewer` valida netlist de exemplo para `examples/` e README.
4. Gates: `pytest && pytest --cov=minispice --cov-fail-under=75 && ruff check src tests scripts tools`.

### Fluxo de biblioteca de células (digital)
1. `cell-designer` entrega a célula (netlist + registro em `tools/gates.py` + truth table + README).
2. `cell-verification` aprova: truth table completa, corners, Monte Carlo yield, noise margins.
3. `cell-char` mede: VTC, CGG, tpHL/tpLH vs CL, energia, tabela por corner.
4. Gates (os mesmos).

Nota DFT: scan/ATPG não existem neste fluxo (sem HDL/síntese aqui); o equivalente didático é o
yield de fabricação do Monte Carlo do `cell-verification` + análise de testabilidade sob demanda.
