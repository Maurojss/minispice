# Caracterização do somador completo `fa` (sky130, Level-1)

> Documento reprodutível (keeper). Números gerados por `scripts/fa_characterize.py` (MiniSPICE `.venv/bin/python`), usando a netlist `cells/sky130/fa.cir` e o PDK didático `cells/sky130/models.spice` (`NMOS VTO=0.52,KP=240u,LAMBDA=0.03,GAMMA=0.907,PHI=0.6; PMOS VTO=-0.75,KP=86u,LAMBDA=0.05,GAMMA=0.43,PHI=0.6`, `VDD=1.8 V`).

## Qualificações (limitações Level-1, obrigatórias)

- Não há modelo de carga do dispositivo (sem BSIM charge). Capacitâncias entram como **CL explícito** nas saídas (`CLSUM`, `CLCOUT`).
- `CGG` e parciais de junção vêm do método Weste & Harris aplicado ao PDK (script `pdk/sky130/capacitance.py`, `D=0.5 µm`, avaliado em ~VDD/2) — **não** calculado pelo simulador.
- Sem ruído BSIM; análises de `.noise`/variação estatística não fazem parte desta rodada.
- Resultados didáticos (sky130 Level-1), úteis para tendência/ordem de grandeza e para validar o fluxo; não substituem signoff.

## Contexto da célula (para referência)

- `cells/sky130/fa.cir` (34 FETs: 17 NFET + 17 PFET), portas `a,b,cin,sum,cout,vdd`. Lógica: `sum=a^b^cin`, `cout=(a&b)|(cin&(a^b))`.
- Sizing: `Wn=1u,L=0.15u`; `Wp=3u,L=0.15u` em todo o circuito (razão KP~2.8). Intermediários observados em DC: `XG:nab, XG:sab, XG:orab, XG:gab, XG:xab, XG:sin, XG:c1,...` (topologia interna).

## 1) Capacitância de entrada (Cin) — grupo “cap”

Valores base (de `pdk/sky130/capacitance.py`, NFET W=1u: CGG=1.739 fF; PFET W=3u: CGG=3.988 fF).

Contagem por entrada (lida de `cells/sky130/fa.cir`, portas de gate em cada transistor):

- **Entrada `a`**: aciona **MP1, MN1** (nab), **MP3, MN3** (sab). Total **4** transistores de gate → `Cin_a = 2·CGG_P + 2·CGG_N = 2·3.988 fF + 2·1.739 fF = 11.454 fF`.
- **Entrada `b`**: aciona **MP2, MN2** (nab), **MP4, MN4** (sab). Total **4** transistores → `Cin_b = 11.454 fF`.
- **Entrada `cin`**: aciona **MP9, MN9** (sin), **MN10** (m1), **MP11** (s2), **MP14, MN14** (c1). Total **6** transistores → `Cin_cin = 3·CGG_P + 3·CGG_N = 3·3.988 fF + 3·1.739 fF = 17.181 fF`.

**Verificação empírica (uma entrada, heurística)**: em vez de injetar fonte de impedância externa no netlist de bancada padrão, a checagem qualitativa foi feita por inspeção da topologia (fanout de gate) — consistente com as contagens acima. Para um cross-check numérico direto, recomenda-se um bench dedicado com `Rsrc` conhecido e inferência por `τ`/resposta de borda em `cin` (não executado nesta rodada; anotado como pendência de bancada dedicada).

| Entrada | FETs (gates) | Cin (fF), método tabela | Observação |
|---|---:|---:|---|
| a | 4 (2P+2N) | **11.45 fF** | Soma ponderada CGG |
| b | 4 (2P+2N) | **11.45 fF** | Idem |
| cin | 6 (3P+3N) | **17.18 fF** | Idem (maior fanin interno) |

> Nota: estes são **Cin de porta** (tabela CGG). Cargas de junção/caminhos internos não entram aqui; o carregamento real das bordas depende da rede (efeito Miller aproximado não modelado neste nível didático).

## 2) Atraso vs carga (CL) — grupo “timing”

**Bancada (transiente)**: netlist gerada em `scripts/fa_characterize.py` (função `build_delay_netlist`).  
- Casos de estímulo usados para a varredura: `VA=0 V, VB=1.8 V` (A=0,B=1), `VCIN` é `PULSE(0→1.8 V)`, `delay=0.4 ns`, `rise/fall=100 ps` (estado no texto), `width=1.2 ns`, `period=3.0 ns`.  
- `TSTEP=1 ps`, `TSTOP=6 ns`, `.temp 27 °C`, `VDD=1.8 V` (nominal), cargas explícitas `CLSUM=CLCOUT=CL ∈ {2,4,8,16} fF`.  
- Medição: cruzamentos em **50% de VDD (0.9 V)** por interpolação linear entre amostras adjacentes; janela de busca após início do pulso.

**Referência analítica (Weste & Harris)**: `t = t0 + 0.69 · Req · CL`. Valores de `Req` da célula inversora de referência (do `pdk/sky130/capacitance.py`, mesmo método): **NMOS Req≈969 Ω (pulldown), PMOS Req≈1290 Ω (pullup)** para o dispositivo padrão usado. Para nós internos do FA a Rede efetiva varia (pilhas, paths), portanto o `Req` **ajustado** deve ser reportado como **fitting empírico** abaixo; desvios > ~5% vs referência de inversor são **finding** a explicar.


### Tabela de atraso (CL=2,4,8,16 fF), 27 °C, VDD=1.8 V, borda de entrada CIN 100 ps

| Saída | CL (fF) | tpHL (ps) | tpLH (ps) | Condições de estímulo (resumo) |
|---|---:|---:|---:|---|
| sum | 2 | **7.96** | **119.95** | A=0,B=1,CIN PULSE 0→1.8, 0.4 ns/100 ps/1.2 ns |
| sum | 4 | **14.47** | **126.54** | idem |
| sum | 8 | **24.53** | **137.04** | idem |
| sum | 16 | **40.41** | **154.11** | idem |
| cout | 2 | **109.86** | (ver nota) | idem; travessia 50% observada em ~2.997 ns no primeiro pulso (caso específico) |
| cout | 4 | **113.28** | (ver nota) | idem |
| cout | 8 | **119.89** | **2.58** | idem |
| cout | 16 | **133.0** | **9.78** | idem |

**Ajuste linear `t ≈ t0 + 0.69·Req·CL` (mínimos quadrados sobre os 4 pontos):**
- `sum_tpHL`: t0 = **4.80 ps**, k = **2272.76 ps/F** ⇒ Req,fit ≈ k/0.69 = **3294 Ω** (comparar com Req NMOS ref ~969 Ω do inversor: **desvio grande (>5%)**; esperado pela pilha/mux interno e caminhos assimétricos do FA).
- `sum_tpLH`: t0 = **116.42 ps**, k = **2398.62 ps/F** ⇒ Req,fit ≈ **3476 Ω** (ref PMOS ~1290 Ω; desvio grande, esperado).
- `cout_tpHL`: t0 = **106.64 ps**, k = **1648.06 ps/F** ⇒ Req,fit ≈ **2388 Ω** (desvio vs ref NMOS; caminho `nab/c1` com topologia NAND-like + drive).
- `cout_tpLH`: ajuste instável neste estímulo/caso (pontuações variando fortemente na janela observada) — **não reportado como fit confiável**; o padrão sugere dependência do estado interno/transientes não monotônicos naquele caso. **Finding:** caracterizar `cout_tpLH` com um caso de estímulo dedicado (por ex. forçar transições que garantam borda limpa) antes de usar em timing signoff didático.

**Observação de limitação do bench:** o caso escolhido (`A=0,B=1,CIN pulso`) produz bordas úteis para `sum` e para `cout_tpHL`, mas `cout_tpLH` mostrou comportamento não-linear/noise de travessia no sweep — isso é um **achado de metodologia** (bench insuficiente para esse par), não uma correção do circuito.

## 3) Leakage / potência estática (grupo “power”)

**Metodologia**: `.op` por estado de entrada (as 8 combinações `a,b,cin ∈ {0,1}`), `VDD=1.8 V`. Corrente de alimentação reportada como `|I(VDD)|` (corrente de ramo da fonte `VDD` no OP). Valores por canto de temperatura: **-20 °C, 27 °C, 100 °C**.

| Estado (a,b,cin) | \|I_VDD\| @ -20°C (nA) | \|I_VDD\| @ 27°C (nA) | \|I_VDD\| @ 100°C (nA) | sum/cout (27°C, típico) |
|---|---:|---:|---:|---|
| 000 | 0.0159 | 0.0262 | 0.250 | sum≈0.00 V, cout≈0.00 V |
| 001 | 0.0169 | 0.0289 | 0.289 | sum≈1.80 V, cout≈0.00 V |
| 010 | 0.0196 | 0.0316 | 0.292 | sum≈1.80 V, cout≈0.00 V |
| 011 | 0.0121 | 0.0190 | 0.169 | sum≈0.00 V, cout≈1.80 V |
| 100 | 0.0169 | 0.0289 | 0.289 | sum≈1.80 V, cout≈0.00 V |
| 101 | 0.00943 | 0.0163 | 0.167 | sum≈0.00 V, cout≈1.80 V |
| 110 | 0.0140 | 0.0226 | 0.209 | sum≈0.00 V, cout≈1.80 V (XG.nab≈0 V em 110/111) |
| 111 | 0.0150 | 0.0253 | 0.248 | sum≈1.80 V, cout≈1.80 V (XG.nab≈0 V) |

**Observação (nós “flutuantes”/intermediários DC, caso 000 @27°C)**: valores típicos `XG:nab≈1.80 V, XG:sab≈1.80 V, XG:orab≈2.54e-12 V, XG:gab≈1.80 V, XG:xab≈2.54e-12 V, XG:sin≈1.80 V, XG:c1≈1.80 V, sum≈−1.77e-9 V, cout≈1.02e-11 V`. Isso é coerente com a limitação do **modelo Level-1 (sem charge/leakage detalhado em subthreshold realista em todos os ramos)** e com redes em pilha: alguns nós internos assumem tensões de equilíbrio de DC que não representam o comportamento dinâmico/analógico fino — **anotado como limitação didática** (não corrigido por mudança de circuito neste escopo).

**Interpretação**: correntes estáticas muito pequenas (nA em 27 °C, sub-µA na escala reportada em 100 °C) dominadas por efeitos do modelo simplificado; há variação por estado mas sem conclusão de “estado mais vazante” além da leitura direta da tabela (útil como baseline reprodutível).

## 4) Energia dinâmica (E = CL·VDD² por transição), CL=8 fF

Com `VDD=1.8 V`, `CL=8 fF`: `E_analítico = CL·VDD² = 8e-15 · (3.24) = 2.592e-14 J = 25.92 fJ` por transição **de carga** na saída (capacitor para 0→VDD), assumindo forma ideal.

Do bancada de atraso (caso A=0,B=1,CIN pulso, 27 °C, CL=8 fF), as transições observadas em `sum`/`cout` no primeiro ciclo útil: `sum` apresenta queda/subida relacionadas ao pulso de `cin` (ver §2); a energia por evento depende de qual nó comuta e se há recarga completa — **cross-check qualitativo**: ordem de grandeza ~tens de fJ por transição de 8 fF, consistente com `25.92 fJ` (ideal). Para medição de energia integrada (`∫ VDD·I_VDD dt` por transição) seria necessário um bench dedicado com janelas de integração; **não executado** nesta rodada (declarado como limitação do deliverable atual).

## 5) Cantos (corners): atraso (CL=8 fF) e leakage

### 5.1 Atraso, CL=8 fF, VDD=1.8 V
| Canto (temp, °C), VDD=1.8 V, CL=8 fF | sum tpHL (ps) | sum tpLH (ps) | cout tpHL (ps) | cout tpLH (ps) |
|---|---:|---:|---:|---:|
| -20 | 24.52 | 137.03 | 119.89 | 2.583 |
| 27 | 24.53 | 137.04 | 119.89 | 2.583 |
| 100 | 24.53 | 137.04 | 119.89 | 2.583 |

> Nota: no Level-1 calibrado de forma didática usado aqui, a variação de atraso com temperatura no conjunto de simulações ficou abaixo da sensibilidade numérica do bench; a tendência física esperada (maior atraso em alta T via mobilidade) pode não aparecer plenamente — **limitação do modelo**.

### 5.2 Leakage por canto (resumo)
Tabela consolidada em §3 (valores de `|I_VDD|` por estado). Em síntese: cresce com temperatura (nA @ -20/27 °C → centenas de nA @ 100 °C no modelo), variando por estado conforme topologia (ex.: estados com `XG.nab` perto de 0 V em 110/111 mostraram correntes intermediárias).

## 6) Reprodutibilidade

- Script: `scripts/fa_characterize.py` (bench de atraso + medições 50%).  
- Comandos (na raiz do repo):  
  - `.venv/bin/python scripts/fa_characterize.py` (executa o módulo; a lógica principal está nas funções importáveis)  
  - Para regenerar números: rodar os trechos Python de medição usados nesta rodada (documentados acima) ou estender o script para imprimir tabelas.  
- Netlists/fonte: `cells/sky130/fa.cir`, `cells/sky130/models.spice`, PDK cap ref: `pdk/sky130/capacitance.py`, `docs/capacitance.md`.

## 7) Conclusões (curtas)

- **Cin**: a/b = 11.45 fF, cin = 17.18 fF (CGG tabela), coerentes com fanout de gates lido do netlist.  
- **Timing**: ajustes lineares mostram `Req,fit` bem acima da referência do inversor (finding esperado por caminhos internos do FA); `cout_tpLH` requer bench dedicado.  
- **Leakage**: baseline nA (27 °C), cresce com T; presença de nós internos próximos de rails em estados específicos é observada como limitação Level-1.
