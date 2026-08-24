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
