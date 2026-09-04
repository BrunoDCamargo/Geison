# Geison

Pipeline para desenho e avaliação in silico de ensaios qPCR/RT-qPCR.

O desenvolvimento ativo acontece na branch `develop`.

## Workbench guiado no Google Colab

Para pesquisadores e avaliadores, use [`notebooks/geison_guided_colab.ipynb`](notebooks/geison_guided_colab.ipynb). O fluxo guiado mantém a lógica científica no Geison e organiza a análise em etapas visíveis:

```text
Target conservation -> Target vs non-target contrast -> Assay design -> Specificity
```

O notebook oferece um **Demo (synthetic)** determinístico e um modo **Project** para arquivos FASTA fornecidos pelo pesquisador. O painel precisa ser revisado e aprovado antes da execução científica. Os resultados são lidos dos artefatos TSV, JSON e HTML produzidos pelo CLI, e o estado final vem de `run_manifest.json`.

O notebook [`notebooks/geison_colab.ipynb`](notebooks/geison_colab.ipynb) continua disponível para validação operacional de baixo nível. Consulte [`docs/guided-colab.md`](docs/guided-colab.md) para o fluxo guiado e [`docs/colab.md`](docs/colab.md) para a operação no Colab.

## Reprodutibilidade e diagnóstico

Antes de executar um pipeline, o ambiente pode ser inspecionado sem configuração e sem acesso à rede:

```bash
qpcr-pipeline doctor
```

Para validar a configuração e visualizar as decisões `RUN`, `REUSE` e `FORCED` sem executar etapas científicas ou gravar artefatos:

```bash
qpcr-pipeline run config.yaml --dry-run
qpcr-pipeline run config.yaml --dry-run --outdir run
```

Uma execução com `--outdir` registra `run_manifest.json` e `run.log.jsonl`, além dos resumos já publicados pelo pipeline. O estado final da run é `COMPLETED`, `PARTIAL` ou `FAILED`. Uma run `PARTIAL` nunca pode publicar `IN SILICO PASS`; evidência incompleta força revisão/incompletude no ranking. Em uma retomada, o mesmo `run_id` é preservado e uma nova tentativa é acrescentada ao histórico.

O `doctor` reporta BLAST+ como `NOT_USED` e não obrigatório. O `--dry-run` não consulta o NCBI, não cria ou modifica o diretório de saída e não executa CD-HIT, MAFFT ou Primer3. Os resultados continuam sendo evidência in silico e não substituem validação experimental do ensaio.

Detalhes sobre manifesto, log sanitizado, proveniência, estados finais e retomada estão em [`docs/reproducible-runs.md`](docs/reproducible-runs.md).

## Checkpoints e retomada

Toda execução com `--outdir` grava checkpoints internos por etapa em
`<outdir>/.checkpoints/`. Uma execução normal continua recalculando todas as
etapas; os checkpoints apenas deixam essa saída pronta para uma retomada futura.

```bash
qpcr-pipeline run config.yaml --outdir run1
qpcr-pipeline run config.yaml --outdir run1 --resume
qpcr-pipeline run config.yaml --outdir run1 --from-step conservation
qpcr-pipeline run config.yaml --outdir run1 --resume --force-step specificity
```

`--resume` valida o manifest, o estado tipado e os SHA-256 dos artefatos de cada
etapa. Checkpoints válidos são reutilizados; uma etapa ausente, incompleta,
corrompida ou incompatível com os inputs, parâmetros ou versões atuais é
recalculada junto com seus dependentes. A invalidação segue o grafo do pipeline,
portanto mudar somente um parâmetro de especificidade não refaz alinhamento ou
conservação, enquanto mudar o clustering invalida toda a cadeia que depende dele.

`--from-step` é estrito: a etapa escolhida e seus dependentes são recalculados,
mas todos os checkpoints necessários fora desse subgrafo precisam continuar
válidos. Se algum pré-requisito estiver inválido, o comando falha antes de iniciar
nova computação científica e informa a etapa bloqueante.

`--force-step` só pode ser usado junto com `--resume`. Ele força a etapa escolhida
e todos os seus dependentes, mas mantém ramos independentes reutilizáveis quando
seus checkpoints continuam válidos. Por exemplo, forçar inclusividade não obriga
a refazer especificidade; o ranking é refeito porque depende das duas evidências.

Os fingerprints incluem apenas os parâmetros relevantes da etapa, identidades dos
inputs, resultados das dependências e a versão do Geison. Quando usados, CD-HIT,
MAFFT e Primer3 também entram com sua identidade de versão somente na etapa que os
invoca. Alteração ou remoção de um output declarado, de `state.json` ou do manifest
torna o checkpoint não reutilizável.

`.checkpoints/` é infraestrutura interna e local ao próprio `outdir`; não existe
cache global. Copiar um diretório de saída completo preserva seu estado local de
retomada, desde que os arquivos continuem íntegros. `run_summary.json` registra
`stage_actions` com `RUN`, `REUSE` ou `FORCED` para mostrar exatamente o que foi
recalculado ou reaproveitado.

Sem `--outdir`, o comando continua apenas carregando e validando a configuração,
e controles de retomada não são aceitos.

## Panel approval workflow

O painel de alvo, grupos obrigatórios e não alvos é uma entrada científica
explícita. O fluxo de aprovação é:

1. Configure `panel.proposal` e execute com `--outdir`.
2. O Geison retorna `ACTION_REQUIRED` e grava `panel_proposal.yaml`.
3. Revise ou edite a proposta como uma entrada científica.
4. Congele a proposta revisada:

   ```bash
   qpcr-pipeline panel approve panel_proposal.yaml --output approved_panel.json
   ```

5. Substitua `panel.proposal` na configuração:

   ```yaml
   panel:
     frozen_manifest: approved_panel.json
   ```

6. Retome a execução:

   ```bash
   qpcr-pipeline run config.yaml --outdir run1 --resume
   ```

O manifesto congelado faz parte da proveniência científica e deve ser versionado
ou arquivado junto com as entradas da execução. A construção automática do painel
não faz parte deste subprojeto.

## Clustering do Discovery Set

O clustering é opcional e vem desabilitado por padrão. Quando habilitado, a
configuração pode incluir identidade, threads e memória:

```yaml
clustering:
  enabled: true
  identity: 0.95
  threads: 4
  memory_mb: 2048
```

`identity` aceita valores de `0.80` a `1.0`. O comprimento de palavra e o
comprimento mínimo efetivo são derivados desse valor; com clustering habilitado,
qualquer sequência aprovada menor que o comprimento de palavra é rejeitada antes da
execução do `cd-hit-est`.

Execuções com clustering habilitado exigem `cd-hit-est` no `PATH`. O Evaluation
Set continua sendo toda a população aprovada pelo QC; o Discovery Set é apenas o
subconjunto de representantes usado para etapas de descoberta.

Os artefatos gerados no diretório de saída são:

- `discovery_set.fasta`: sequências representantes do Discovery Set.
- `clustering_report.json`: configuração, contagens, membros dos clusters e
  rastreabilidade entre Evaluation Set e Discovery Set.
- `clustering/cd-hit-est.clstr`: saída bruta do CD-HIT, presente apenas quando o
  clustering está habilitado; uma nova execução desabilitada remove esse artefato
  obsoleto sem remover outros arquivos do diretório `clustering/`.

## Alinhamento do Discovery Set

O alinhamento é opcional e vem desabilitado por padrão. Quando habilitado, ele
alinha somente os representantes do Discovery Set, após o clustering, e pode ser
configurado assim:

```yaml
alignment:
  enabled: true
  threads: 4
  reference_id: seq-3
```

`reference_id` é opcional. Quando informado, deve identificar uma sequência do
Discovery Set, e essa sequência é usada como referência. Sem `reference_id`, a
referência é escolhida automaticamente pela menor fração de bases ambíguas,
depois pela maior sequência e, em caso de empate, pela ordem do Discovery Set.

Execuções com alinhamento habilitado e duas ou mais sequências no Discovery Set
exigem `mafft` no `PATH`. O MAFFT recebe `--adjustdirectionaccurately`; quando
ele identifica uma sequência em orientação reversa complementar, o resultado
normalizado é publicado com o ID original e a orientação registrada no relatório.

Os artefatos de alinhamento são:

- `alignment/alignment_report.json`: configuração, referência, orientações,
  contagens e caminhos dos artefatos.
- `alignment/discovery_alignment.fasta`: sequências alinhadas, na ordem do
  Discovery Set.
- `alignment/coordinate_map.tsv`: mapa entre a coluna do alinhamento e a
  coordenada/base da referência.

`qc_report.json` também inclui `alignment` com status, ID de referência e modo de
seleção (`explicit` ou `automatic`) para rastreabilidade de alto nível. Ambas as
colunas de coordenadas são 1-based; gaps na referência deixam vazios os campos de
referência.

## Conservação genômica

A análise de conservação é opcional, vem desabilitada por padrão e depende de um
alinhamento habilitado. A janela e o passo podem ser configurados assim:

```yaml
alignment:
  enabled: true
  reference_id: seq-3

conservation:
  enabled: true
  window_size: 100
  step_size: 10
```

`window_size` e `step_size` aceitam inteiros de 1 a 1.000.000, com o passo menor
ou igual à janela. A validação ocorre antes da criação do diretório de saída.

As métricas são calculadas em todas as colunas do alinhamento do Discovery Set:

- profundidade é o número de sequências sem gap na coluna;
- cobertura é a profundidade dividida pelo número de sequências;
- frequência de gap é o número de gaps dividido pelo número de sequências;
- códigos IUPAC ambíguos distribuem um voto fracionário igualmente entre suas
  bases A, C, G e T compatíveis; gaps não participam dessas frequências;
- conservação é a frequência do alelo majoritário;
- entropia é a entropia de Shannon, em bits, das frequências A/C/G/T;
- o consenso majoritário prefere a base canônica da referência em empates e,
  depois, a ordem A, C, G, T; o segundo consenso preserva o suporte por IUPAC.

As colunas de inserção, nas quais a referência tem gap, permanecem nas métricas
por posição com a coordenada de referência vazia. Elas não entram nos consensos
nem nas janelas. As janelas usam coordenadas 1-based da referência e registram a
média e o mínimo de conservação, além das médias de cobertura, gaps e entropia.
Uma referência menor que a janela produz uma única janela parcial; referências
maiores recebem janelas completas e uma janela terminal ancorada quando necessário.

Quando a referência vem de GenBank, features locais válidas são convertidas para
intervalos 1-based inclusivos e exibidas como anotações. Features `source`, partes
externas e intervalos fora da referência são ignorados de forma rastreável.

Os artefatos publicados são:

- `conservation/conservation_report.json`: configuração, definições, contagens e
  caminhos dos artefatos;
- `conservation/position_metrics.tsv`: métricas por coluna do alinhamento;
- `conservation/window_metrics.tsv`: métricas agregadas por janela;
- `conservation/consensus_major.fasta`: consenso majoritário sem inserções da
  referência;
- `conservation/consensus_iupac.fasta`: consenso IUPAC sem inserções da referência;
- `report.html`: relatório Canvas autocontido, sem CDN ou recursos de rede.

O relatório permite zoom pela roda do mouse, pan por arraste, restauração da visão
completa, detalhes por hover e zoom ao clicar nas janelas mais conservadas. Esta
etapa apenas calcula e visualiza os picos. Se o ranking final estiver habilitado
na mesma execução, ele assume depois a propriedade de `report.html` e substitui
este relatório pelo relatório consolidado dos assays.

Com conservação desabilitada, apenas
`conservation/conservation_report.json` é publicado com status `SKIPPED`; dados
científicos e `report.html` não são gerados. `qc_report.json` sempre inclui o
status, a referência e as contagens de posições e janelas da etapa.

## Desenho de ensaios com Primer3

O desenho de primers e sondas é opcional, vem desabilitado por padrão e depende
de alinhamento e conservação habilitados. A forma YAML completa da seção, com
seus valores padrão, é:

```yaml
primer_design:
  enabled: false
  max_candidate_regions: 10
  assays_per_region: 5
  candidate_region_length: 300
  max_region_overlap_fraction: 0.5
  min_mean_conservation: 0.90
  min_minimum_conservation: 0.70
  min_mean_coverage: 0.90
  max_mean_gap_frequency: 0.05
  max_mean_entropy_bits: 0.50
  min_usable_fraction: 0.80
  product_size_min: 70
  product_size_max: 200
  primer:
    min_size: 18
    opt_size: 20
    max_size: 25
    min_tm: 58.0
    opt_tm: 60.0
    max_tm: 62.0
    min_gc_percent: 40.0
    max_gc_percent: 60.0
  probe:
    min_size: 18
    opt_size: 25
    max_size: 30
    min_tm: 68.0
    opt_tm: 70.0
    max_tm: 72.0
    min_gc_percent: 30.0
    max_gc_percent: 80.0
```

Cada janela de conservação elegível é expandida para
`candidate_region_length`; nas extremidades, o intervalo é deslocado para caber
na referência, e uma referência menor usa seu comprimento total. Intervalos
iguais são deduplicados. Os candidatos são ordenados por maior conservação
média, maior conservação mínima, maior cobertura média, menor entropia média,
menor frequência média de gaps, maior comprimento utilizável e, por fim,
coordenadas iniciais e finais menores. Depois de aceitar um candidato, outro só
é aceito quando o comprimento de sua interseção dividido pelo comprimento do
menor dos dois intervalos é menor ou igual a
`max_region_overlap_fraction`. Os IDs `region-001`, `region-002`, etc. seguem
essa ordem.

Quando existe ao menos uma região candidata, a execução exige o executável
`primer3_core` no `PATH`. O consenso majoritário completo é enviado em cada
`SEQUENCE_TEMPLATE`, enquanto `SEQUENCE_INCLUDED_REGION` limita o desenho ao
candidato. Os artefatos ficam em `primer_design/`:

- `primer_design_report.json`: configuração efetiva, contagens, diagnósticos do
  Primer3, ensaios e caminhos dos artefatos;
- `candidate_regions.tsv`: regiões candidatas ordenadas e suas métricas;
- `assays.tsv`: pares completos de primer forward, sonda interna e primer
  reverse;
- `primer3_input.txt` e `primer3_output.txt`: Boulder-IO enviado e recebido para
  auditoria, presentes somente quando o Primer3 é executado.

Todas as coordenadas publicadas nos TSV e JSON são 1-based e inclusivas, tanto
para candidatos quanto para primers e sondas. O tamanho do produto corresponde
a `reverse_reference_end - forward_reference_start + 1`. O Boulder-IO bruto é
preservado sem conversão nos dois artefatos de auditoria.

Com `primer_design.enabled: false`, somente
`primer_design/primer_design_report.json` é publicado para a etapa, com status
`SKIPPED`; o runner não é invocado. Com a etapa habilitada mas sem regiões
elegíveis, o status é `COMPLETE`, `candidate_regions.tsv` e `assays.tsv` contêm
somente os cabeçalhos, e o Primer3 não é executado. Se o Primer3 não retornar
pares completos, `assays.tsv` fica apenas com o cabeçalho e a entrada, a saída e
os diagnósticos continuam preservados. `qc_report.json` inclui o status, a
referência e as contagens de candidatos e ensaios.

Esta etapa produz candidatos de ensaio auditáveis.

## Inclusividade e propostas IUPAC

A avaliação de inclusividade é opcional, vem desabilitada por padrão e depende
do desenho de ensaios habilitado. A configuração efetiva completa é:

```yaml
inclusivity:
  enabled: false
  search_flank: 250
  max_hits_per_oligo: 20
  max_primer_mismatches: 2
  max_probe_mismatches: 1
  reject_primer_3_prime_mismatch: true
  primer_3_prime_bases: 5
  max_primer_degeneracy: 16
  max_probe_degeneracy: 4
  allow_primer_3_prime_degeneracy: false
  max_amplicon_size_delta: 20
```

Cada ensaio do Primer3 é avaliado contra todas as sequências aprovadas do
Evaluation Set, inclusive membros que não foram escolhidos como representantes
do Discovery Set. Cada sequência é pesquisada tanto na orientação fornecida como
na sua reversa complementar. Oligos IUPAC usam uma comparação conservadora: uma
base-alvo ambígua só é considerada coberta quando todas as suas bases possíveis
estão no suporte do símbolo do oligo.

As posições de mismatch são 1-based na orientação de síntese 5-prime para
3-prime do oligo, inclusive para o primer reverse. Por padrão, qualquer mismatch
nos cinco nucleotídeos da extremidade 3-prime de um primer torna esse hit
incompatível. A compatibilidade completa exige a geometria estrita
forward < probe < reverse e um tamanho de amplicon dentro de
`max_amplicon_size_delta` do tamanho projetado.

Quando uma variação observada melhora a cobertura exata e respeita os limites de
degenerescência, o pipeline publica uma proposta IUPAC limitada. O ensaio
original permanece imutável e é avaliado lado a lado com a proposta; nenhuma
proposta substitui silenciosamente o oligo original. Os artefatos ficam em
`inclusivity/`:

- `oligo_matches.tsv`: hits locais e detalhes de mismatch por oligo.
- `assay_inclusivity.tsv`: geometria e compatibilidade original/proposta por
  ensaio e sequência.
- `oligo_variations.tsv`: variações posicionais observadas no Evaluation Set.
- `degeneracy_proposals.tsv`: sequências originais, propostas, limites e motivos.
- `inclusivity_report.json`: relatório normalizado que referencia os quatro TSVs.

As propostas IUPAC são candidatas computacionais para auditoria. Elas não
substituem validação experimental, nem estimam sozinhas consequências
termodinâmicas ou risco biológico.

A issue #8 avalia inclusividade contra todo o Evaluation Set. A issue #9 avalia
especificidade contra conjuntos off-target; propostas IUPAC nunca substituem
silenciosamente os oligos originais. A decisão final de risco e a interface de
usuário pertencem a etapas posteriores.

## Especificidade contra off-targets

A análise de especificidade é opcional, vem desabilitada por padrão e depende do
desenho de ensaios habilitado. Ela não depende de a etapa de inclusividade estar
habilitada. Quando `specificity.enabled: true`, deve existir pelo menos um dataset
off-target configurado, cada um com nome único e exatamente uma fonte.

A configuração efetiva é:

```yaml
off_targets:
  - name: human
    fasta: data/human_subset.fasta
  - name: near_neighbors
    frozen_dataset: runs/ncbi_near_neighbors

specificity:
  enabled: false
  max_hits_per_oligo_per_dataset: 20
  max_primer_mismatches: 2
  max_probe_mismatches: 1
  reject_primer_3_prime_mismatch: true
  primer_3_prime_bases: 5
  max_amplicon_size: 1000
```

Os datasets podem ser FASTA locais ou datasets NCBI já congelados. A etapa de
especificidade nunca consulta o NCBI e não usa BLAST no MVP. Para datasets
congelados, o manifest existente é validado e referenciado para preservar query,
accessions e versões materializadas. Para FASTA local, o relatório registra
caminho, SHA-256 e IDs das sequências.

Forward, probe e reverse são pesquisados exaustivamente na sequência fornecida e
na reversa complementar. A comparação IUPAC segue a mesma semântica conservadora
da inclusividade: uma base-alvo ambígua só é coberta quando todo o suporte desse
símbolo está contido no suporte do símbolo do oligo. Oligos degenerados são
comparados diretamente, sem expansão combinatória obrigatória.

A etapa diferencia três situações:

- um hit isolado de oligo, sem geometria F/R válida;
- `primer_amplicon_plausible`, quando forward e reverse compatíveis estão voltados
  um para o outro e delimitam um intervalo de até `max_amplicon_size`;
- `detectable_off_target`, quando esse amplicon plausível também contém uma probe
  compatível entre os primers.

`max_hits_per_oligo_per_dataset` limita somente os hits individuais publicados em
`off_target_hits.tsv`. Todos os hits compatíveis são avaliados antes desse limite
para formar a geometria; portanto o truncamento do TSV não pode transformar um
off-target real em resultado aparentemente seguro. Os amplicons preservam as
coordenadas dos primers e probes que sustentam a classificação mesmo quando um
desses hits não aparece entre os hits individuais retidos.

Os artefatos ficam em `specificity/`:

- `off_target_hits.tsv`: hits individuais retidos, com dataset, assay, sequência,
  papel do oligo, orientação, coordenadas e mismatches;
- `plausible_amplicons.tsv`: todas as geometrias F/R plausíveis e o estado de
  detecção pela probe;
- `specificity_report.json`: configuração efetiva, proveniência dos datasets,
  contagens, truncamentos e caminhos dos artefatos.

Com a etapa desabilitada, somente `specificity_report.json` é publicado com
status `SKIPPED` e nenhum dataset off-target é lido. `qc_report.json` inclui um
resumo com quantidade de datasets, sequências, assays, hits retidos, amplicons
plausíveis e off-targets detectáveis.

A busca exaustiva em Python é destinada a conjuntos off-target pequenos ou
moderados, como near-neighbors e coleções de referência congeladas e curadas. Ela
não foi projetada para varrer bancos genômicos muito grandes. Se essa necessidade
surgir, um backend indexado, como BLAST, poderá ser adicionado sem mudar o contrato
científico de hits, geometria e detectabilidade definido nesta etapa.

## Classificação e ranking final dos assays

O ranking final é opcional, vem desabilitado por padrão e depende do desenho de
ensaios habilitado. Inclusividade e especificidade podem estar desabilitadas, mas
nesse caso a ausência de evidência impede `IN SILICO PASS` e o assay fica pelo
menos em `REVIEW` com reason codes explícitos.

A configuração padrão é:

```yaml
ranking:
  enabled: false
  min_inclusivity_for_pass: 1.0
  min_inclusivity_before_high_risk: 0.90
  weights:
    inclusivity: 0.35
    specificity: 0.25
    conservation: 0.20
    primer3_quality: 0.10
    robustness: 0.10
```

A classificação acontece antes do score. Com os valores padrão, inclusividade
original de 100% não rebaixa a classe; valores de 90% até menos de 100% geram
`REVIEW`; abaixo de 90% gera `HIGH_RISK`. Na especificidade, qualquer
`detectable_off_target` gera `HIGH_RISK`; um amplicon F/R plausível sem probe
compatível gera `REVIEW`; hits isolados são apenas advisory e reduzem o componente
de especificidade.

O score é absoluto, determinístico e vai de 0 a 100. Ele é decomposto em cinco
componentes nomeados: inclusividade, especificidade, conservação, qualidade
Primer3 e robustez. Os pesos padrão são, respectivamente, 35%, 25%, 20%, 10% e
10%. A classe é sempre a primeira chave de ordenação, portanto nenhum score alto
permite que `HIGH_RISK` ultrapasse `REVIEW` ou `IN SILICO PASS`.

Se qualquer componente necessário estiver indisponível, `score_status` passa a
`INCOMPLETE` e `final_score` fica vazio, mesmo que o peso configurado daquele
componente seja zero. Evidência ausente não é transformada artificialmente em
score zero.

Propostas IUPAC aceitas ou rejeitadas continuam sendo evidência contextual. O
assay original do Primer3 é o único candidato ranqueado, e uma proposta degenerada
não é reavaliada automaticamente contra off-targets nem substitui silenciosamente
F/Probe/R.

Os artefatos são:

- `ranking/assay_ranking.tsv`: uma linha por assay, incluindo classe, score,
  componentes e reason codes;
- `ranking/ranking_report.json`: configuração efetiva, contagens, componentes e
  reasons estruturados com evidências;
- `report.html`: relatório final autocontido com F/Probe/R, métricas da região,
  inclusividade, propostas de degeneração, especificidade, classificação e
  ranking.

Com ranking desabilitado, o estágio não altera `report.html`; portanto o relatório
publicado pela conservação continua disponível. Com ranking habilitado, o ranking
assume a propriedade de `report.html` e substitui o relatório de conservação pelo
relatório final consolidado dos assays.

`IN SILICO PASS` significa apenas que o assay satisfez as regras computacionais e
as evidências disponíveis desta execução. Essa classificação não constitui nem
substitui validação experimental do ensaio.
