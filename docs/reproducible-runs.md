# Execuções reprodutíveis e diagnóstico

A issue #12 adiciona diagnóstico de ambiente, pré-visualização sem efeitos colaterais e registro auditável das execuções do Geison.

## Diagnóstico do ambiente

Use `doctor` antes de uma execução para inspecionar Python, Geison, Git e ferramentas externas sem carregar configuração e sem acessar a rede:

```bash
qpcr-pipeline doctor
```

CD-HIT, MAFFT e Primer3 aparecem como disponíveis ou indisponíveis. A obrigatoriedade efetiva depende das etapas habilitadas na configuração. BLAST+ é reportado como `NOT_USED` e `required: false`; ele não é dependência obrigatória do Geison.

## Dry-run

Para validar a configuração, inspecionar o ambiente e visualizar as decisões `RUN`, `REUSE` e `FORCED` sem executar etapas científicas:

```bash
qpcr-pipeline run config.yaml --dry-run
qpcr-pipeline run config.yaml --dry-run --outdir run
```

O dry-run não cria nem modifica o diretório de saída, não consulta o NCBI e não executa CD-HIT, MAFFT ou Primer3.

## Manifesto e log da execução

Uma execução com `--outdir` publica:

- `run_manifest.json`: identidade da run, configuração efetiva sanitizada, ambiente, política de execução, plano, histórico append-only de tentativas, proveniência das entradas, referência, completude científica e falha final quando aplicável.
- `run.log.jsonl`: eventos estruturados da execução, uma linha JSON por evento, com campos permitidos e sanitização de informações sensíveis.
- `run_summary.json`: resumo final da tentativa concluída.
- `qc_report.json`: resumo das etapas científicas e suas contagens.

Ao iniciar uma nova tentativa no mesmo diretório, `run_summary.json` e `qc_report.json` anteriores são invalidados para evitar que um resultado antigo pareça ser o resultado da tentativa atual. Os checkpoints não são apagados e continuam disponíveis para `--resume`.

## Estados finais

O estado da run é um destes:

- `COMPLETED`: todas as evidências científicas necessárias para o fluxo habilitado estão completas.
- `PARTIAL`: a execução terminou, mas faltam evidências científicas necessárias para uma conclusão completa.
- `FAILED`: a tentativa falhou e a falha foi persistida no manifesto.

Uma run `PARTIAL` nunca pode publicar um ensaio como `IN SILICO PASS`. Evidência incompleta força revisão/incompletude no ranking, preservando `HIGH_RISK` quando já existe evidência de alto risco.

Em uma retomada, o `run_id` permanece o mesmo e uma nova tentativa é acrescentada ao histórico do manifesto. Tentativas anteriores não são reescritas.

## Proveniência

Para entradas locais, o manifesto reutiliza a identidade SHA-256 já calculada pelo checkpoint de entrada e registra formato e contagens de QC. O arquivo não é re-hashado por um caminho paralelo.

Para entradas NCBI, o manifesto projeta somente os campos necessários para auditoria: modo de aquisição, query ou accessions solicitados, versões de accession resolvidas e identidade do dataset. Metadados internos de batches, payloads brutos, e-mail NCBI e chave de API não são copiados para o manifesto de run.

## Limites da avaliação in silico

Os resultados do Geison são evidência computacional para seleção e revisão de candidatos. `IN SILICO PASS` não substitui validação experimental, avaliação analítica nem os controles aplicáveis ao uso pretendido do ensaio.
