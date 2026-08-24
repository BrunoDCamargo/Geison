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
etapa apenas calcula e visualiza os picos. A selecao de regioes candidatas para
desenho de ensaios pertence a issue #7.

Com conservacao desabilitada, apenas
`conservation/conservation_report.json` e publicado com status `SKIPPED`; dados
cientificos e `report.html` nao sao gerados. `qc_report.json` sempre inclui o
status, a referencia e as contagens de posicoes e janelas da etapa.
