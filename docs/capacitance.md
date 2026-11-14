# Capacitância de parasita dos transistores — método Weste & Harris no sky130

Exemplo completo de **caracterização capacitiva** dos dispositivos usados nas
nossas células compostas, seguindo a metodologia do livro **CMOS VLSI
Design (Weste & Harris, cap. 4 — "Circuit Characterization")** e usando
parâmetros **reais** dos modelos BSIM4 do sky130 PDK vendidos em
`pdk/sky130/models/`.

Dispositivos (os mesmos das células):

- NFET: `W = 1.0 µm`, `L = 0.15 µm`
- PFET: `W = 3.0 µm`, `L = 0.15 µm`
- `VDD = 1.8 V`, corner `tt`

Todos os números abaixo saem de `python pdk/sky130/capacitance.py`.

## 1. Capacitor de porta (CGG)

O livro modela a porta como o capacitor Cox em série com a camada de óxido:

```
CGG = Cox · W · L  +  2 · CGDO · W          (canal + overlaps fonte/dreno)
Cox = ε_ox / tox,   ε_ox = 3.9 · ε0
```

Do PM3 do sky130 (`cgso = cgdo`):

| parâmetro | NFET | PFET | unidade |
|---|---|---|---|
| `toxm` | 4.148 | 4.230 | nm |
| Cox | 8.325 | 8.163 | mF/m² |
| `cgso/cgdo` | 2.449 | 0.525 | nF/m |

Resultado:

| dispositivo | canal | overlap | **CGG total** | fF/µm |
|---|---|---|---|---|
| NFET 1.0 µm | 1.25 fF | 0.49 fF | **1.74 fF** | 1.74 |
| PFET 3.0 µm | 3.67 fF | 0.31 fF | **3.99 fF** | 1.33 |

Comparação com a regra de bolso do livro (≈ 2 fF/µm em processos antigos):
sky130 dá **1.7 fF/µm** para o NFET — coerente para um processo moderno
de óxido fino. (A capacidade fringa/porta-bulk `cgbo = 0.1 nF/m` é desprezível:
~0.0001 fF nos nossos tamanhos.)

## 2. Capacitores de junção (Csb/Cdb)

Modelo do Harris com coeficientes de graduação:

```
Cj(V) = CJ · AD · (1 + V/PB)^(−MJ)          (fundamental, area)
      + CJSW  · PD · (1 + V/PBSW)^(−MJSW)   (sidewall lateral)
      + CJSWG · PG · (1 + V/PBSWG)^(−MJSWG) (gate-edge sidewall)
```

Parâmetros BSIM4 equivalentes no sky130 (`cjs`, `cjsws`, `cjswgs`, `mjs`,
`mjsws`, `mjswgs`, `pbs`, `pbsws`, `pbswgs`):

| | CJ [mF/m²] | MJ | PB | CJSW [nF/m] | MJSW | PBSW |
|---|---|---|---|---|---|---|
| NFET | 1.340 | 0.44 | 0.729 | 0.0367 | 0.0009 | 0.200 |
| PFET | 0.738 | 0.346 | 0.659 | 0.0989 | 0.298 | 0.742 |

Geometria do strip de difusão (hipótese documentada, altura `D = 0.5 µm`):

```
AD = W · D        PD = 2·D + W     PG = W
```

> Ajuste `--D` ao seu layout; a capacitância cresce linearmente com AD e PD.

Avaliadas em `V = VDD/2` (método do livro: valor eficaz no meio da varredura):

| dispositivo | Cj(0) | **Cj(VDD/2)** |
|---|---|---|
| NFET | 0.98 fF | **0.68 fF** |
| PFET | 2.22 fF | **1.59 fF** |

## 3. Carga efetiva do inversor (FO1)

```
C_L = self-loading + próximas portas = (0.68 + 1.59) + (1.74 + 3.99) = 8.00 fF
```

## 4. Resistência efetiva e atraso estimado

Framework do livro: `t = 0.69 · Req · C_L`, com

```
Req = (1/0.69) · ∫_{VDD/2}^{VDD}  dV / Id(Vds=V)
```

usando as próprias equações MOS do simulador (nível 1 com λ, calibradas ao
bin do PDK — ver `pdk/sky130/summary.md`):

| | Req | previsão 0.69·Req·CL |
|---|---|---|
| pulldown (NMOS) | 969 Ω | **tPHL = 5.35 ps** |
| pull-up (PMOS) | 1290 Ω | **tPLH = 7.12 ps** |

## 5. Verificação em simulação (MiniSPICE)

`examples/cmos_inverter_cap.cir` aplica o pulso e carrega a saída com
`CL = 8 fF` explícito:

| | 0.69·Req·CL (livro) | MiniSPICE (.tran) | erro |
|---|---|---|---|
| tPHL | 5.35 ps | 5.61 ps | 4.9% |
| tPLH | 7.12 ps | 7.41 ps | 4.1% |

O erro < 5% valida o fluxo completo: **PDK real → Level-1 → capacitância
→ previsão de atraso → simulação**. As diferenças restantes vêm do modelo
(sem carga de fiação, Miller e não linearidades BSIM4).

## Como reproduzir / adaptar

```bash
python pdk/sky130/capacitance.py            # tabela completa + previsões
python pdk/sky130/capacitance.py --D 0.3    # outro layout de difusão
minispice examples/cmos_inverter_cap.cir    # verificação por simulação
```

Para outros dispositivos, mude `W`, `L` e `KP/VTO/LAMBDA` na chamada `load()`
do script — os parâmetros BSIM4 (tox, cg, cj...) são extraídos direto dos
arquivos reais em `models/`.
