# Geison

Pipeline para desenho e avaliacao in silico de ensaios qPCR/RT-qPCR.

O desenvolvimento ativo acontece na branch `develop`.

## Clustering da Discovery Set

O clustering e opcional e vem desabilitado por padrao. Quando habilitado, a
configuracao pode incluir identidade, threads e memoria:

```yaml
clustering:
  enabled: true
  identity: 0.95
  threads: 4
  memory_mb: 2048
```

`identity` aceita valores de `0.80` a `1.0`. O comprimento de palavra e o
comprimento minimo efetivo sao derivados desse valor; com clustering habilitado,
qualquer sequencia aprovada menor que o comprimento de palavra e rejeitada antes da
execucao do `cd-hit-est`.

Execucoes com clustering habilitado exigem `cd-hit-est` no `PATH`. A Evaluation
Set continua sendo toda a populacao aprovada pelo QC; a Discovery Set e apenas o
subconjunto de representantes usado para etapas de descoberta.

Os artefatos gerados no diretorio de saida sao:

- `discovery_set.fasta`: sequencias representantes da Discovery Set.
- `clustering_report.json`: configuracao, contagens, membros dos clusters e
  rastreabilidade entre Evaluation Set e Discovery Set.
- `clustering/cd-hit-est.clstr`: saida bruta do CD-HIT, presente apenas quando o
  clustering esta habilitado; uma nova execucao desabilitada remove esse artefato
  obsoleto sem remover outros arquivos do diretorio `clustering/`.

## Alinhamento da Discovery Set

O alinhamento e opcional e vem desabilitado por padrao. Quando habilitado, ele
alinha somente os representantes da Discovery Set, apos o clustering, e pode ser
configurado assim:

```yaml
alignment:
  enabled: true
  threads: 4
  reference_id: seq-3
```

`reference_id` e opcional. Quando informado, deve identificar uma sequencia da
Discovery Set e essa sequencia e usada como referencia. Sem `reference_id`, a
referencia e escolhida automaticamente pela menor fracao de bases ambiguas,
depois pela maior sequencia e, em caso de empate, pela ordem da Discovery Set.

Execucoes com alinhamento habilitado e duas ou mais sequencias na Discovery Set
exigem `mafft` no `PATH`. O MAFFT recebe `--adjustdirectionaccurately`; quando
ele identifica uma sequencia em orientacao reversa complementar, o resultado
normalizado e publicado com o ID original e a orientacao registrada no relatorio.

Os artefatos de alinhamento sao:

- `alignment/alignment_report.json`: configuracao, referencia, orientacoes,
  contagens e caminhos dos artefatos.
- `alignment/discovery_alignment.fasta`: sequencias alinhadas, na ordem da
  Discovery Set.
- `alignment/coordinate_map.tsv`: mapa entre a coluna do alinhamento e a
  coordenada/base da referencia.

`qc_report.json` tambem inclui `alignment` com status, ID de referencia e modo de
selecao (`explicit` ou `automatic`) para rastreabilidade de alto nivel. Both
coordinate columns are 1-based while reference gaps have blank reference fields.

## Conservacao genomica

A analise de conservacao e opcional, vem desabilitada por padrao e depende de um
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
ou igual a janela. A validacao ocorre antes da criacao do diretorio de saida.

As metricas sao calculadas em todas as colunas do alinhamento da Discovery Set:

- profundidade e o numero de sequencias sem gap na coluna;
- cobertura e a profundidade dividida pelo numero de sequencias;
- frequencia de gap e o numero de gaps dividido pelo numero de sequencias;
- codigos IUPAC ambiguos distribuem um voto fracionario igualmente entre suas
  bases A, C, G e T compativeis; gaps nao participam dessas frequencias;
- conservacao e a frequencia do alelo majoritario;
- entropia e a entropia de Shannon, em bits, das frequencias A/C/G/T;
- o consenso majoritario prefere a base canonica da referencia em empates e,
  depois, a ordem A, C, G, T; o segundo consenso preserva o suporte por IUPAC.

As colunas de insercao, nas quais a referencia tem gap, permanecem nas metricas
por posicao com a coordenada de referencia vazia. Elas nao entram nos consensos
nem nas janelas. As janelas usam coordenadas 1-based da referencia e registram a
media e o minimo de conservacao, alem das medias de cobertura, gaps e entropia.
Uma referencia menor que a janela produz uma unica janela parcial; referencias
maiores recebem janelas completas e uma janela terminal ancorada quando necessario.

Quando a referencia vem de GenBank, features locais validas sao convertidas para
intervalos 1-based inclusivos e exibidas como anotacoes. Features `source`, partes
externas e intervalos fora da referencia sao ignorados de forma rastreavel.

Os artefatos publicados sao:

- `conservation/conservation_report.json`: configuracao, definicoes, contagens e
  caminhos dos artefatos;
- `conservation/position_metrics.tsv`: metricas por coluna do alinhamento;
- `conservation/window_metrics.tsv`: metricas agregadas por janela;
- `conservation/consensus_major.fasta`: consenso majoritario sem insercoes da
  referencia;
- `conservation/consensus_iupac.fasta`: consenso IUPAC sem insercoes da referencia;
- `report.html`: relatorio Canvas autocontido, sem CDN ou recursos de rede.

O relatorio permite zoom pela roda do mouse, pan por arraste, restauracao da visao
completa, detalhes por hover e zoom ao clicar nas janelas mais conservadas. Esta
etapa apenas calcula e visualiza os picos.

Com conservacao desabilitada, apenas
`conservation/conservation_report.json` e publicado com status `SKIPPED`; dados
cientificos e `report.html` nao sao gerados. `qc_report.json` sempre inclui o
status, a referencia e as contagens de posicoes e janelas da etapa.

## Desenho de ensaios com Primer3

O desenho de primers e sondas e opcional, vem desabilitado por padrao e depende
de alinhamento e conservacao habilitados. A forma YAML completa da secao, com
seus valores padrao, e:

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

Cada janela de conservacao elegivel e expandida para
`candidate_region_length`; nas extremidades, o intervalo e deslocado para caber
na referencia, e uma referencia menor usa seu comprimento total. Intervalos
iguais sao deduplicados. Os candidatos sao ordenados por maior conservacao
media, maior conservacao minima, maior cobertura media, menor entropia media,
menor frequencia media de gaps, maior comprimento utilizavel e, por fim,
coordenadas iniciais e finais menores. Depois de aceitar um candidato, outro so
e aceito quando o comprimento de sua intersecao dividido pelo comprimento do
menor dos dois intervalos e menor ou igual a
`max_region_overlap_fraction`. Os IDs `region-001`, `region-002`, etc. seguem
essa ordem.

Quando existe ao menos uma regiao candidata, a execucao exige o executavel
`primer3_core` no `PATH`. O consenso majoritario completo e enviado em cada
`SEQUENCE_TEMPLATE`, enquanto `SEQUENCE_INCLUDED_REGION` limita o desenho ao
candidato. Os artefatos ficam em `primer_design/`:

- `primer_design_report.json`: configuracao efetiva, contagens, diagnosticos do
  Primer3, ensaios e caminhos dos artefatos;
- `candidate_regions.tsv`: regioes candidatas ordenadas e suas metricas;
- `assays.tsv`: pares completos de primer forward, sonda interna e primer
  reverse;
- `primer3_input.txt` e `primer3_output.txt`: Boulder-IO enviado e recebido para
  auditoria, presentes somente quando o Primer3 e executado.

Todas as coordenadas publicadas nos TSV e JSON sao 1-based e inclusivas, tanto
para candidatos quanto para primers e sondas. O tamanho do produto corresponde
a `reverse_reference_end - forward_reference_start + 1`. O Boulder-IO bruto e
preservado sem conversao nos dois artefatos de auditoria.

Com `primer_design.enabled: false`, somente
`primer_design/primer_design_report.json` e publicado para a etapa, com status
`SKIPPED`; o runner nao e invocado. Com a etapa habilitada mas sem regioes
elegiveis, o status e `COMPLETE`, `candidate_regions.tsv` e `assays.tsv` contem
somente os cabecalhos, e o Primer3 nao e executado. Se o Primer3 nao retornar
pares completos, `assays.tsv` fica apenas com o cabecalho e a entrada, a saida e
os diagnosticos continuam preservados. `qc_report.json` inclui o status, a
referencia e as contagens de candidatos e ensaios.

Esta etapa produz candidatos de ensaio auditaveis.

## Inclusividade e propostas IUPAC

A avaliacao de inclusividade e opcional, vem desabilitada por padrao e depende
do desenho de ensaios habilitado. A configuracao efetiva completa e:

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

Cada ensaio do Primer3 e avaliado contra todas as sequencias aprovadas do
Evaluation Set, inclusive membros que nao foram escolhidos como representantes
da Discovery Set. Cada sequencia e pesquisada tanto na orientacao fornecida como
na sua reversa complementar. Oligos IUPAC usam uma comparacao conservadora: uma
base alvo ambigua so e considerada coberta quando todas as suas bases possiveis
estao no suporte do simbolo do oligo.

As posicoes de mismatch sao 1-based na orientacao de sintese 5-prime para
3-prime do oligo, inclusive para o primer reverse. Por padrao, qualquer mismatch
nos cinco nucleotideos da extremidade 3-prime de um primer torna esse hit
incompativel. A compatibilidade completa exige a geometria estrita
forward < probe < reverse e um tamanho de amplicon dentro de
`max_amplicon_size_delta` do tamanho projetado.

Quando uma variacao observada melhora a cobertura exata e respeita os limites de
degenerescencia, o pipeline publica uma proposta IUPAC limitada. O ensaio
original permanece imutavel e e avaliado lado a lado com a proposta; nenhuma
proposta substitui silenciosamente o oligo original. Os artefatos ficam em
`inclusivity/`:

- `oligo_matches.tsv`: hits locais e detalhes de mismatch por oligo.
- `assay_inclusivity.tsv`: geometria e compatibilidade original/proposta por
  ensaio e sequencia.
- `oligo_variations.tsv`: variacoes posicionais observadas na Evaluation Set.
- `degeneracy_proposals.tsv`: sequencias originais, propostas, limites e motivos.
- `inclusivity_report.json`: relatorio normalizado que referencia os quatro TSVs.

As propostas IUPAC sao candidatas computacionais para auditoria. Elas nao
substituem validacao experimental, nem estimam sozinhas consequencias
termodinamicas ou risco biologico.

A issue #8 avalia inclusividade contra todo o Evaluation Set. A issue #9 avalia
especificidade contra conjuntos off-target; propostas IUPAC nunca substituem
silenciosamente os oligos originais. A decisao final de risco e a interface de
usuario pertencem a etapas posteriores.
