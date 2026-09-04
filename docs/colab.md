# Google Colab

Para pesquisadores e avaliadores, o fluxo recomendado está em [`notebooks/geison_guided_colab.ipynb`](../notebooks/geison_guided_colab.ipynb), com instruções detalhadas em [`docs/guided-colab.md`](guided-colab.md).

Ele apresenta a sequência científica como:

```text
Target conservation -> Target vs non-target contrast -> Assay design -> Specificity
```

O notebook guiado oferece modo **Demo (synthetic)** e modo **Project**, mantém o gate humano `APROVAR`, mostra os estados de `run_manifest.json` e renderiza apenas artefatos publicados pelo Geison. A lógica científica continua no pacote e no CLI.

O notebook anterior, [`notebooks/geison_colab.ipynb`](../notebooks/geison_colab.ipynb), permanece disponível para validação operacional de baixo nível e troubleshooting do ambiente.

Ambos preparam o ambiente, escrevem a configuração YAML e chamam o CLI `qpcr-pipeline`. O fluxo inclui o gate de painel antes da execução científica: `panel.proposal` -> `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` -> revisão humana -> `panel approve` -> `frozen_manifest` -> `--resume`.

## Primeira execução

Em uma sessão nova, o notebook clona a branch `main` e instala o pacote em modo editável:

```bash
git clone --branch main https://github.com/BrunoDCamargo/Geison.git
python -m pip install -e /content/Geison
```

Se `/content/Geison` já existir, ele atualiza o checkout:

```bash
git -C /content/Geison checkout main
git -C /content/Geison pull --ff-only origin main
python -m pip install -e /content/Geison
```

Executado de dentro do repositório, o comando equivalente de atualização é `git pull --ff-only origin main`.

O notebook imprime o SHA em teste com `git rev-parse HEAD`, instala CD-HIT, MAFFT e Primer3 e executa:

```bash
qpcr-pipeline doctor
```

## Identificação para aquisição NCBI

Aquisição NCBI ao vivo exige `NCBI_EMAIL`. `NCBI_API_KEY` é opcional. O notebook mantém os dois valores somente no ambiente da sessão; não os grava no YAML, manifestos ou relatórios.

## Gate de painel

A primeira configuração contém `panel.proposal`. O `--dry-run` deve anunciar que a aprovação é necessária sem criar o diretório de saída:

```bash
qpcr-pipeline run /content/geison_run/config.yaml --dry-run --outdir /content/geison_run/output
```

A primeira execução real para com código de processo `3`, status `ACTION_REQUIRED` e código `PANEL_APPROVAL_REQUIRED`. Ela grava `/content/geison_run/output/panel_proposal.yaml` e não deve criar o checkpoint de `input`.

O notebook exibe o `panel_proposal.yaml` para revisão. A aprovação exige uma confirmação humana explícita (`APROVAR`) antes de executar:

```bash
qpcr-pipeline panel approve \
  /content/geison_run/output/panel_proposal.yaml \
  --output /content/geison_run/approved_panel.json
```

O manifesto aprovado deve conter `status: APPROVED`, `approved_by_user: true` e `proposal_sha256` com prefixo `sha256:`.

Depois da aprovação, a configuração passa a usar:

```yaml
panel:
  frozen_manifest: /content/geison_run/approved_panel.json
```

A mesma run é retomada no mesmo diretório:

```bash
qpcr-pipeline run /content/geison_run/config-approved.yaml \
  --outdir /content/geison_run/output \
  --resume
```

Após o `--resume`, o notebook confirma:

- `panel/approved_panel.json`;
- `.checkpoints/panel/manifest.json`;
- `.checkpoints/input/manifest.json`;
- `run_manifest.json` com `panel_provenance`.

## Relatórios

O relatório específico de contraste fica em:

```text
/content/geison_run/output/contrastive_conservation/report.html
```

Quando a configuração publica o relatório final consolidado, ele fica em:

```text
/content/geison_run/output/report.html
```

No notebook guiado, as visualizações de conservação, contraste, desenho, cobertura, especificidade e ranking são produzidas a partir dos TSV/JSON/HTML publicados. `PARTIAL`, `FAILED` e `ACTION_REQUIRED` não são apresentados como conclusão bem-sucedida.

## Persistência no Colab

O armazenamento em `/content` é temporário. Para retomar em outra sessão, preserve `config-approved.yaml`, `approved_panel.json` e o diretório completo `/content/geison_run/output`. Copiar apenas os arquivos finais não preserva os checkpoints exigidos por `--resume`.
