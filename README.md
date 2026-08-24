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
