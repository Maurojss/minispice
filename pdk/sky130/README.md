# sky130 PDK (vendido neste repositório)

Esta pasta contém cópias literais dos **modelos BSIM4 reais** (corner `tt`,
PM3) dos transistores `nfet_01v8` e `pfet_01v8` do **SkyWater sky130 PDK**,
junto com o script que extrai os parâmetros usados para calibrar os modelos
Level-1 em `cells/sky130/models.spice`.

```
pdk/sky130/
├── models/nfet_tt.pm3.spice   # BSIM4 (level=54, version 4.5), 180 bins, corner tt
├── models/pfet_tt.pm3.spice   # idem, PFET
├── extract.py                # extrai os parâmetros do bin de interesse
└── summary.md                # sumário dos parâmetros e o mapeamento Level-1
```

- **Fonte:** `google/skywater-pdk-libs-sky130_fd_pr` (GitHub), arquivos
  `cells/<dispositivo>/<dispositivo>__tt.pm3.spice`. Licença Apache-2.0.
- **O que são PM3:** os modelos de 3ᵃ geração da SkyWater — BSIM4.5 completo,
  com 180 "bins" (um `.model` por faixa de L×W) por dispositivo.
- **Por que não usamos estes arquivos na simulação:** eles exigem um SPICE
  com BSIM4 (`level=54`) — ex.: ngspice compilado com o flag de BSIM4. O
  MiniSPICE implementa MOS **Level-1** (Shichman–Hodges); em vez de fingir
  que os modelos BSIM4 rodam aqui, extraímos os parâmetros reais do bin
  escolhido e usamos equivalentes Level-1 calibrados por eles. O resultado é
  qualitativamente correto (ponto de comutação ~VDD/2, razões de drive
  coerentes), não é um sign-off.

Para conferir os números extraídos do bin usado pelas células:

```bash
python pdk/sky130/extract.py
```

Os binários escolhidos e o mapeamento completo estão em [summary.md](summary.md).
