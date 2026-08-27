# Issue #9: especificidade contra off-targets

Status: Approved design
Date: 2026-08-27
Issue: #9 `Executar especificidade contra off-targets`
Branch: `feature/issue-9-specificity`

## Objetivo

Adicionar ao Geison uma etapa de especificidade in silico que avalie assays qPCR contra múltiplos datasets off-target de forma offline, determinística e reproduzível. A etapa deve distinguir um hit isolado de oligo, um amplicon plausível formado por primers e um off-target potencialmente detectável pela sonda.

O MVP não usa BLAST. A busca é interna e exaustiva sobre FASTA local ou datasets NCBI previamente congelados. Um backend escalável para bancos genômicos grandes fica fora do escopo e poderá ser introduzido futuramente sem alterar o contrato científico da etapa.

## Decisões principais

1. Criar um módulo `qpcr_pipeline/specificity.py` independente.
2. Reutilizar apenas primitives públicas e estáveis, especialmente `iupac.py` e os modelos de assay do desenho de primers.
3. Não reutilizar funções privadas de `inclusivity.py` e não generalizar a etapa #8 durante esta issue.
4. Executar especificidade depois da inclusividade no `run_pipeline()`.
5. Manter `specificity.enabled: false` por padrão para preservar compatibilidade com configurações existentes.
6. Aceitar somente datasets off-target offline: FASTA local ou dataset NCBI congelado.
7. Tratar `primer_amplicon_plausible` e `detectable_off_target` como estados diferentes.

## Arquitetura

Fluxo de alto nível:

```text
Primer design
    |
Inclusivity
    |
Specificity
    |
    +-- carregar N datasets off-target
    |     +-- FASTA local
    |     +-- NCBI congelado
    |
    +-- pesquisar F / Probe / R em duas orientações
    |
    +-- preservar hits individuais
    |
    +-- combinar hits F/R por geometria
    |
    +-- verificar probe dentro do amplicon
    |
    +-- publicar artefatos e resumo
```

O módulo de especificidade terá seus próprios tipos públicos, previstos inicialmente como:

- `OffTargetDataset`
- `OffTargetHit`
- `PlausibleAmplicon`
- `SpecificityResult`
- `SpecificityError`

Os nomes exatos podem ser ajustados durante a implementação desde que os conceitos permaneçam separados e o contrato dos artefatos seja preservado.

## Configuração

Os datasets off-target ficam separados dos parâmetros científicos da etapa.

Exemplo:

```yaml
off_targets:
  - name: human
    fasta: data/human_subset.fasta

  - name: near_neighbors
    frozen_dataset: runs/ncbi_near_neighbors

specificity:
  enabled: true
  max_hits_per_oligo_per_dataset: 20
  max_primer_mismatches: 2
  max_probe_mismatches: 1
  reject_primer_3_prime_mismatch: true
  primer_3_prime_bases: 5
  max_amplicon_size: 1000
```

### Regras de configuração

- `off_targets` pode conter zero ou mais datasets.
- Cada dataset precisa de `name` não vazio e exatamente uma fonte: `fasta` ou `frozen_dataset`.
- Nomes de datasets devem ser únicos.
- `specificity.enabled: true` exige que `primer_design.enabled: true`.
- A etapa pode ser habilitada sem `inclusivity.enabled`; a ordem no pipeline continua sendo depois da inclusividade quando ambas estiverem habilitadas.
- `max_hits_per_oligo_per_dataset` deve ser inteiro positivo.
- `max_primer_mismatches` e `max_probe_mismatches` devem ser inteiros não negativos.
- `primer_3_prime_bases` deve ser inteiro positivo.
- `max_amplicon_size` deve ser inteiro positivo e começa em 1000 bp.

## Origem e proveniência dos datasets

### FASTA local

O dataset preserva:

- nome configurado;
- caminho efetivo;
- SHA-256 do arquivo;
- quantidade de sequências;
- IDs das sequências.

### NCBI congelado

A etapa usa somente dados já materializados e validados pelo mecanismo de frozen dataset existente. Não há acesso de rede durante a especificidade.

O relatório referencia o manifest do dataset congelado e preserva a composição materializada, incluindo query, accessions e accession versions quando disponíveis no manifest existente.

## Semântica de matching

A busca é exaustiva e determinística sobre cada sequência off-target e sua reversa complementar.

### IUPAC

A especificidade usa a mesma semântica conservadora da inclusividade:

- um símbolo IUPAC do oligo representa o conjunto de bases que ele aceita;
- uma base-alvo ambígua só é considerada coberta quando todo o seu conjunto de possibilidades está contido no suporte do símbolo do oligo;
- oligos degenerados são comparados diretamente e não precisam ser expandidos combinatoriamente.

A implementação não chama helpers privados de `inclusivity.py`. Um teste de equivalência semântica protege contra divergência entre as duas etapas.

### Mismatches

- posições de mismatch são 1-based na orientação de síntese 5' para 3' do oligo;
- forward, reverse e probe mantêm contagem e posições de mismatch;
- primers respeitam `max_primer_mismatches`;
- probe respeita `max_probe_mismatches`;
- quando `reject_primer_3_prime_mismatch` estiver habilitado, qualquer mismatch nos últimos `primer_3_prime_bases` torna o hit do primer incompatível.

### Ordenação determinística

Antes do truncamento, os hits são ordenados de forma estável por:

`dataset -> assay -> sequence -> role -> orientation -> coordinate -> mismatch_count`

Critérios adicionais de desempate necessários à implementação devem ser explícitos e estáveis.

Somente depois da ordenação é aplicado `max_hits_per_oligo_per_dataset`. Truncamentos são registrados no relatório.

## Geometria do assay

A análise distingue três níveis:

1. **Hit isolado**: um oligo possui um hit compatível, mas não existe geometria F/R válida.
2. **`primer_amplicon_plausible`**: forward e reverse possuem hits compatíveis, estão voltados um para o outro e delimitam um intervalo com tamanho menor ou igual a `max_amplicon_size`.
3. **`detectable_off_target`**: existe um amplicon plausível e uma probe compatível dentro do intervalo delimitado pelos primers.

A presença de probe fora do intervalo não torna o amplicon detectável.

Múltiplas combinações plausíveis podem ser preservadas. A implementação não escolhe silenciosamente apenas uma combinação quando existem várias geometrias válidas.

## Artefatos

A etapa publica:

```text
specificity/
├── off_target_hits.tsv
├── plausible_amplicons.tsv
└── specificity_report.json
```

### `off_target_hits.tsv`

Cada linha representa um hit individual e inclui, no mínimo:

- dataset;
- assay_id;
- sequence_id;
- role;
- orientation;
- coordenadas;
- mismatch_count;
- mismatch_positions;
- compatibilidade;
- informação de truncamento quando aplicável ao conjunto de hits.

### `plausible_amplicons.tsv`

Cada linha representa uma combinação F/R plausível e inclui, no mínimo:

- dataset;
- assay_id;
- sequence_id;
- orientation;
- coordenadas do amplicon;
- tamanho;
- referência aos hits forward e reverse;
- presença e compatibilidade da probe dentro do intervalo;
- `primer_amplicon_plausible`;
- `detectable_off_target`.

### `specificity_report.json`

Consolida:

- status da etapa;
- configuração efetiva;
- proveniência de cada dataset;
- contagens de sequências, hits e amplicons;
- contagens de assays com risco;
- registros de truncamento;
- caminhos dos artefatos.

`qc_report.json` recebe somente um resumo da etapa, seguindo o padrão existente do pipeline.

## Comportamento de erro e estados

### `SKIPPED`

Quando `specificity.enabled: false`:

- nenhum dataset é lido;
- não há busca;
- `specificity_report.json` é publicado com status `SKIPPED`;
- TSVs científicos não precisam ser publicados.

### `COMPLETE`

A etapa termina como `COMPLETE` quando a avaliação termina corretamente, inclusive quando:

- não existem assays;
- um dataset válido contém zero sequências;
- nenhum hit é encontrado.

Nesses casos, artefatos tabulares publicados devem permanecer válidos, com cabeçalhos quando aplicável.

### Falhas explícitas

A etapa falha com erro contextualizado quando houver:

- dataset ausente ou ilegível;
- FASTA inválido;
- frozen dataset inválido;
- símbolo IUPAC inválido;
- configuração inválida;
- inconsistência interna que impeça garantir resultado determinístico.

Erros de IUPAC devem informar, quando disponível, assay, dataset, sequence e role.

## Limites de escala

A busca interna do MVP é intencionalmente simples e determinística. Ela é adequada para datasets off-target pequenos ou moderados usados em desenvolvimento, near-neighbors e conjuntos de referência congelados.

Ela não é indicada para varrer bancos genômicos muito grandes. BLAST ou outro backend indexado será considerado somente quando houver necessidade medida. Um backend futuro deve consumir o mesmo contrato de datasets e produzir semanticamente os mesmos tipos de hits e amplicons, para não alterar a interpretação científica da etapa.

## Testes

A implementação segue TDD. Os testes previstos são:

```text
tests/
├── test_specificity_config.py
├── test_specificity_matching.py
├── test_specificity_geometry.py
├── test_specificity_artifacts.py
└── test_pipeline_specificity.py
```

Cobertura obrigatória:

- ausência de hit;
- hit isolado de forward, reverse e probe;
- F/R em orientação inválida;
- F/R corretos, mas acima de `max_amplicon_size`;
- probe fora do amplicon;
- amplicon plausível sem probe compatível;
- amplicon plausível com probe detectável;
- sequência em orientação reversa complementar;
- mismatch permitido;
- mismatch proibido no 3' de primer;
- IUPAC degenerado;
- base-alvo IUPAC ambígua;
- múltiplos hits com ordenação estável;
- truncamento de hits;
- FASTA válido;
- FASTA vazio;
- FASTA inválido;
- frozen dataset válido;
- frozen dataset inválido;
- especificidade desabilitada;
- nenhum assay desenhado;
- geração dos três artefatos;
- integração `primer_design -> inclusivity -> specificity`;
- equivalência semântica básica de IUPAC entre inclusividade e especificidade.

Os testes de #9 não dependem de rede, BLAST, MAFFT, CD-HIT ou Primer3, exceto testes de integração já existentes que exercitem essas etapas por seus runners/fixtures conforme o padrão do projeto.

## CircleCI

Os testes da #9 entram na suíte normal de `pytest` executada em `develop`.

A branch `feature/issue-9-specificity` não consome executor com a configuração atual, porque o CircleCI executa automaticamente apenas em `develop` e `main`.

Fluxo de integração:

```text
feature/issue-9-specificity
    -> revisão
    -> merge em develop
    -> uma execução CircleCI
    -> PASS: fechar #9
    -> FAIL: corrigir e repetir
```

## Critérios de conclusão da issue #9

A issue só pode ser fechada como `completed` quando:

1. todos os acceptance criteria atualizados da #9 estiverem implementados;
2. os testes específicos da etapa estiverem verdes;
3. a suíte completa de `develop` estiver verde no CircleCI;
4. os artefatos e a limitação de escala estiverem documentados no README;
5. a etapa não fizer acesso de rede;
6. a revisão do diff não encontrar pendência crítica ou importante.

O teste manual do usuário é uma validação adicional de uso e não bloqueia o fechamento técnico da #9. Problemas encontrados nessa validação posterior devem virar bug ou regressão rastreável.

## Fora do escopo

- BLAST+;
- `blastn-short`;
- construção de bases BLAST;
- acesso NCBI ao vivo durante a especificidade;
- indexação para bancos genômicos grandes;
- refatoração geral da inclusividade;
- classificação e ranking final dos assays, que pertencem à #10.
