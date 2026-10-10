# Artigo × código original: divergências encontradas na migração

Este documento registra as diferenças entre o texto publicado de

> de Souza, J.N.M., de Medeiros, J.L., Costa, A.L.H., Nunes, G.C. (2010). *Modeling,
> simulation and optimization of continuous gas lift systems for deepwater offshore
> petroleum production*. Journal of Petroleum Science and Engineering,
> doi:10.1016/j.petrol.2010.03.028 (`docs/paper/PETROL_1857.pdf`)

e o código MATLAB `bbnet10` que gerou os seus números (`legacy/matlab/bbnet10`). A
rastreabilidade entre os arquivos `.mat` e as tabelas/figuras do artigo está em
`legacy/matlab/RESUMO_CASOS_E_ARQUIVOS.txt` e foi confirmada numericamente: os valores
gravados nos `.mat` coincidem com todos os dígitos impressos nas Tabelas 13 a 19.

O pacote `gaslift` **reproduz o código por padrão** (é o que permite validar a migração
contra os resultados originais) e expõe parâmetros para as formas descritas no artigo.
O notebook `06_artigo_vs_codigo.ipynb` quantifica o efeito de cada item.

## 1. Geometria dos poços (Tabelas 2 e 7)

| Grandeza | Artigo | Código (`case_03.m`, `case_04.m`) | Parâmetro em `gaslift.cases` |
|---|---|---|---|
| Diâmetro da coluna de produção | 5 in | 4 in | `wellbore_diameter_in` |
| Comprimento do anular (topo → ponto de injeção) | 1500 m | 1710 m | `annulus_length` |
| Trecho inferior da coluna (injeção → fundo) | 300 m | 338 m (2048 − 1710) | `wellbore_length` |
| Rugosidade das linhas de gás | 100 µm (Tab. 7) / 0 (Tab. 2) | 0 | `gas_line_roughness` |
| Perfil das linhas 20 e 26 (Tabela 10) | 0, 1005, 2505, **4005** m | 0, 1005, 2505, **3005** m | — |

A Tabela 16 usa $L_{G3} = L_{G4} = 4520$ m = 3005 + 1515 m, consistente com o código.

## 2. Orientação geométrica do anular

Em `creating_wells` o caminho do anular é escrito de baixo para cima
(`pipe_z = [z_wellhead − L_annulus, z_wellhead]`), enquanto a aresta do digrafo vai do
topo (nó n+3) ao ponto de injeção (nó n+1). O integrador trata o anular como escoamento
**ascendente**: o termo hidrostático tem o sinal trocado em relação ao gás que desce. A
Tabela 2 do artigo lista −90°, que é a intenção física.

Efeito: o lado do gás não realimenta o lado da produção (as especificações são vazões),
logo pressão de fundo, produção máxima e Tabelas 14–15 não mudam; a pressão de descarga
do compressor fica ≈ 2 × (coluna de gás de 1710 m) ≈ 45 bar maior, afetando potência,
custo de gás combustível e custos de capital. `WellDesign(annulus_orientation="physical")`
corrige; `Network.validate()` acusa a inconsistência na rede padrão.

## 3. Diâmetro equivalente do anular

O código usa $D_{eq} = 4(A_o + A_i)/(P_o + P_i) = (D_o^2 + D_i^2)/(D_o + D_i)$
(8,29 in para 10 in × 4 in). O diâmetro hidráulico de um anular é $D_o − D_i$ (6 in).
`WellDesign(annulus_diameter_model="hydraulic")` corrige.

## 4. Pressão de descarga do compressor ($P_G$)

O índice `ix_Pgaslift` aponta para o nó n+4 (interconexão das linhas de gás na base do
riser, −1495 m) e não para o nó n+5 (entrada de gás na plataforma, +20 m) citado no texto
("pressure at the gas inlet (node 10)"). As Figs. 7–8 e as pressões de descarga das
Tabelas 13 e 18 são, portanto, pressões na base do riser; a pressão na plataforma é menor
pela coluna hidrostática de gás da linha vertical de 1515 m (≈ 17–20 bar).
`CostModel(discharge_pressure_node="gas_inlet")` corrige.

## 5. Custo das linhas de gás (Eq. 33, McCoy 2008)

O código multiplica o diâmetro pelo expoente (`D_nps.*a7`) em vez de elevá-lo
(`D_nps.^a7`). A Eq. 33 impressa está correta; os números da Tabela 18 (CP = 1,05 e
0,98 MUSD) só se reproduzem com a forma do código (a forma correta dá ≈ 1,12 MUSD).
`CostModel(pipeline_cost_formulation="mccoy")` usa a Eq. 33.

## 6. Parâmetros do compressor e correlações de custo

* A Tabela 16 lista coeficiente politrópico 1,25, mas a potência das Tabelas 13 e 18 só se
  reproduz com $\gamma = 1{,}4$ (`parametersOpt.comp.gama = 1.4`).
  `CompressorParameters(polytropic_exponent=1.25)` usa o valor da tabela.
* As Eqs. 30 e 32 (Couper et al.) são escritas com a potência em kW; o código converte
  para hp (`bhp = 1.341 W_c`) antes de aplicá-las, e são esses os valores das tabelas.
  `CostModel(compressor_cost_power_unit="kW")` usa kW.

## 7. Contagem de equações

Caso 2: o artigo informa 232 equações ("27 + 27 + 27 + 56 + 95"); com N = 28 nós e
S = 27 linhas, a formulação da Seção 3.3 tem $3N + 4S = 192$ equações e $N + S = 55$
especificações. O Caso 1 tem 66 equações, como publicado (embora a decomposição impressa,
"9 + 9 + 9 + 20 + 19", também não corresponda a N = 10 e S = 9).

## 8. Tolerância da integração

Os `.mat` das varreduras do Caso 1 (Figs. 6–8) foram gerados com `OdeRelTol = 1e-3`,
enquanto `case_03.m` e todos os demais resultados usam `1e-2`. O integrador em lote do
pacote replica o controle de passo do `ode45`; com a tolerância certa os vetores de estado
coincidem com o MATLAB em ~1e-10.

## 9. Termos cinéticos do balanço de momento nodal (Eq. 13)

A expressão implementada, $|v|\,A\,q/(2\rho^2)$, não é dimensionalmente uma pressão e vale
~1e-14 (em unidades de 100 bar): é numericamente irrelevante. Foi mantida por fidelidade.

## 10. Vazões "nulas" das fontes

Os nós fonte especificam 1e-6 kg/s para a fase ausente (`gasInletwL`, `downholewG`), e
não zero, para afastar a correlação de Beggs & Brill dos limites monofásicos. O pacote
faz o mesmo (`gaslift.cases.EPSILON_FLOW`).
