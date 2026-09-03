# Google Colab

O fluxo oficial do Geison para Google Colab está em [`notebooks/geison_colab.ipynb`](../notebooks/geison_colab.ipynb).

O notebook é deliberadamente fino: ele prepara o ambiente, escreve a configuração YAML e chama o CLI `qpcr-pipeline`. Nenhum algoritmo científico é mantido no notebook; clustering, alinhamento, conservação, desenho de primers, inclusividade, especificidade e ranking continuam pertencendo ao pacote.

## Primeira execução

Abra o notebook no Colab e execute as células em ordem. Em uma sessão nova ele clona a branch `main` e instala o pacote em modo editável:

```bash
git clone --branch main https://github.com/BrunoDCamargo/Geison.git
python -m pip install -e /content/Geison
```

O notebook instala também CD-HIT, MAFFT e Primer3 pelo sistema e, antes da análise, executa:

```bash
qpcr-pipeline doctor
```

Use o resultado do `doctor` para confirmar as versões e a disponibilidade das dependências antes de iniciar uma run.

## Identificação para aquisição NCBI

Aquisição NCBI ao vivo exige a variável de ambiente `NCBI_EMAIL`. O notebook solicita esse e-mail antes da primeira execução científica e o mantém somente no ambiente da sessão:

```python
os.environ["NCBI_EMAIL"] = ncbi_email
```

`NCBI_API_KEY` é opcional. Quando informada, a chave é lida de forma oculta e também fica apenas no ambiente da sessão. E-mail e API key não são colocados no YAML, manifestos ou relatórios do Geison.

## Atualizar uma sessão existente

Se `/content/Geison` já existir, a célula de preparação atualiza o checkout com fast-forward e reinstala o pacote editável:

```bash
git -C /content/Geison checkout main
git -C /content/Geison pull --ff-only origin main
python -m pip install -e /content/Geison
```

O comando equivalente, executado de dentro do repositório, é:

```bash
git pull --ff-only origin main
```

Assim, atualizar o notebook não exige copiar lógica científica nem editar células de implementação.

## Configuração e execução

A análise é controlada por um arquivo `config.yaml`. O notebook inclui um exemplo NCBI pequeno que pode ser substituído pela configuração do projeto.

Valide primeiro sem executar etapas científicas:

```bash
qpcr-pipeline run /content/geison_run/config.yaml --dry-run --outdir /content/geison_run/output
```

Depois execute a run:

```bash
qpcr-pipeline run /content/geison_run/config.yaml --outdir /content/geison_run/output
```

O diretório de saída contém o `run_manifest.json`, o log estruturado, checkpoints e os artefatos científicos habilitados pela configuração.

## Retomar uma run

Para reutilizar checkpoints válidos após interrupção ou reinício da sessão, restaure o mesmo diretório de saída e execute:

```bash
qpcr-pipeline run /content/geison_run/config.yaml --outdir /content/geison_run/output --resume
```

O `--resume` mantém a identidade da run, registra uma nova tentativa e recalcula apenas etapas cujos checkpoints não sejam reutilizáveis.

## `report.html`

Quando a configuração habilita uma etapa que publica o relatório, o arquivo fica em:

```text
/content/geison_run/output/report.html
```

A última célula do notebook lê esse `report.html` e o exibe diretamente no Colab. O mesmo arquivo pode ser baixado ou copiado para o Google Drive.

## Persistência no Colab

O armazenamento em `/content` é temporário. Para uma run longa ou para retomar em outra sessão, preserve `config.yaml` e o diretório completo de saída fora do runtime temporário, por exemplo no Google Drive. Copiar somente arquivos finais não preserva necessariamente os checkpoints necessários ao `--resume`.
