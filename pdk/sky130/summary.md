# sky130 — resumo dos parâmetros e o mapeamento Level-1

Modelos BSIM4 reais do sky130 PDK (corner `tt`) estão em
`models/*_tt.pm3.spice`. Os valores abaixo foram extraídos com `extract.py`
do bin em que as células deste repositório foram formuladas.

## Bins usados

| dispositivo | bin | L | W |
|---|---|---|---|
| `nfet_01v8` | bin pequeno | 0.15 µm – 0.18 µm | 1.00 µm – 1.26 µm |
| `pfet_01v8` | bin pequeno | 0.15 µm – 0.18 µm | 1.00 µm – 1.26 µm |

Parâmetros BSIM4.5 brutos extraídos desses bins:

| parâmetro | BSIM4 | NFET bin | PFET bin |
|---|---|---|---|
| vth0 | `vth0` | 3.19e-2¹ | −0.7956 |
| u0 | `u0` | 2.8796e-2 | −4.91e-3² |
| k1 | `k1` | 0.9071 | 4.118 |
| ndep | `ndep` | 1.7e17 | 1.7e17 |
| toxm | `toxm` | 4.148 nm | 4.23 nm |

¹ vth0 do bin pequeno entra na fórmula `{vth0 + slope/sqrt(LW)}` com termos
de nível (`lvth0`, `wvth0`) — o valor "terminal" é o 0.52 típico do
dispositivo (`vth0 = 0.5190093` no bin de W grande; é o valor adotado).
² u0 do PFET varia por bin (−4.9e-3 … +1.05e-2; os bins grandes valem
+1.0477e-2) — o sinal/polaridade é tratado pelo modelo; em Level-1 adotamos
|u0| típico ≈ 1.05e-2.

## Mapeamento para Level-1 (Shichman–Hodges)

O MiniSPICE resolve o MOS Level-1: `Id = ½·KP·(W/L)·(Vgs−Vt)²·(1+λVds)` com
efeito de corpo `Vt = VTO + γ(√(φ−Vbs) − √φ)`. Correspondências usadas:

| Level-1 | expressão | NFET | PFET | origem |
|---|---|---|---|---|
| `VTO` | vth0 do bin panorâmico | **0.52 V** | **−0.75 V** | 0.519 (NFET); PFET varia −0.41…−1.06 por bin — adotado valor médio típico |
| `KP` | u0·Cox, Cox = εox/tox | **240 µA/V²** | **86 µA/V²** | Cox,NFET = 8.34 mF/m² (tox 4.148 nm); Cox,PFET = 8.16 mF/m² (tox 4.23 nm) |
| `GAMMA` | k1 | **0.907** | **0.43** | k1 bin (NFET) e valor panorâmico do PFET (0.4345) |
| `PHI` | 2·φF(ndep) | **0.6 V** | **0.6 V** | ndep = 1.7e17 ⇒ φF ≈ 0.30 V |
| `LAMBDA` | — (didático) | 0.03 | 0.05 | valor pedagógico; BSIM4 modela de outra forma |

## Escolha das dimensões das células

A razão de correntes a mesma (W/L):

    KP_n / KP_p = 240µ / 86µ ≈ 2.8  ⇒  Wp ≈ 3·Wn

Adotado em todas as células (`cells/sky130/`):

- NFET: `W = 1 µm`, PFET: `W = 3 µm`;
- `L = 0.15 µm` para ambos (L mínimo real desses dispositivos no PDK; para a
  rede PMOS isso é uma idealização didática, já que séries de 2–4 PFETs em
  VDD real teriam W maior para casar resistência).

Teste com `VDD = 1.8 V` (alimentação nominal do processo 1.8 V do sky130).

## Verificação

Com esses parâmetros, todas as 8 células compostas verificam a tabela-verdade
completa em `.op` (66 combinações: `tests/test_compound_gates.py`) e a VTC do
AOI21 com `b=1, c=0` comuta em **VA ≈ 0.82 V** (próximo de VDD/2), plot em
`docs/aoi21_vtc.png`.
