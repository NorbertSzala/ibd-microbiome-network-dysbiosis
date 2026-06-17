# Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes

The project compares two views of the gut microbiome in inflammatory bowel disease (IBD):

- **taxonomic profiles**: which microbes are present,
- **functional pathway profiles**: which metabolic functions are present.

The main question is whether IBD is better reflected by microbial composition, functional potential, or both.

---

## Data

The workflow uses public metagenomic profiles from [`curatedMetagenomicData`](https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html?utm_source=chatgpt.com).

Default dataset in `config/config.yaml`:

```yaml
dataset:
  study_name: "HMP_2019_ibdmdb"
  body_site: "stool"
  disease_column: "study_condition"
```

The pipeline uses two feature matrices:

| Matrix               | Meaning                                        |
| -------------------- | ---------------------------------------------- |
| `taxa_matrix.csv`    | microbial taxonomic abundance profiles         |
| `pathway_matrix.csv` | microbial metabolic pathway abundance profiles |

The final comparison is simplified to two groups:

- `healthy`
- `IBD`

---

## Workflow overview

The workflow contains seven main steps:

| Step | Rule                     | Main output                                          |
| ---- | ------------------------ | ---------------------------------------------------- |
| 01   | `download_data`          | raw `SummarizedExperiment` objects and metadata      |
| 02   | `preprocess`             | matched and filtered taxa/pathway matrices           |
| 03   | `abundance_distribution` | quality-control abundance plots                      |
| 04   | `diversity_analysis`     | alpha-diversity metrics and tests                    |
| 05   | `ordination`             | PCA, PCoA, PERMANOVA results                         |
| 06   | `differential_abundance` | taxa/pathways associated with IBD status             |
| 07   | `classification`         | random forest and elastic net classification results |



---

## Repository structure

```text
.
├── README.md
├── Snakefile
├── Results
|   ├──figures
├── config/
│   └── config.yaml
├── docs/
│   ├── README.md
│   ├── configuration.md
│   ├── extended_readme.md
│   ├── methods.md
│   ├── project_map.md
│   └── workflow_overview.md
├── envs/
│   ├── r_analysis.yaml
│   ├── r_classification.yaml
│   ├── r_download_data.yaml
│   └── snakemake.yml
├── rules/
│   ├── 01_download_data.smk
│   ├── 02_preprocess.smk
│   ├── 03_abundance_distribution.smk
│   ├── 04_diversity_analysis.smk
│   ├── 05_ordination.smk
│   ├── 06_differential_abundance.smk
│   └── 07_classification.smk
├── schemas/
│   └── config.schema.yaml
└── scripts/
    ├── 01_download_data.R
    ├── 02_preprocess.R
    ├── 03_abundance_distribution.R
    ├── 04_diversity_analysis.R
    ├── 05_ordination.R
    ├── 06_differential_abundance.R
    ├── 07_classification.R
    └── functions/
```

Data and result files are generated during workflow execution and are not meant to be stored in the repository file size limitations.

---

## Quick start

Create and activate the Snakemake environment:

```bash
conda env create -f envs/snakemake.yml
conda activate snakemake
```

Run a dry run first:

```bash
snakemake -n --use-conda
```

Run the workflow:

```bash
snakemake --cores 4 --use-conda
```



The first full run downloads data from `curatedMetagenomicData` and creates processed input matrices before running the analysis steps.

NOTE: At the end of finishing this job, the `curatedMetagenomicData` had some troubles with their servers so the data download and package installation was not possible.

---

## Interpretation

The project compares IBD signal across two feature sets:

- stronger taxa results suggest that microbial composition captures more disease signal,
- stronger pathway results suggest that functional potential captures more disease signal,
- similar results suggest that both views contain related information.

Classification results should be interpreted together with diversity, ordination, and differential abundance results. They should not be treated as a clinical diagnostic model.

---

## Documentation

More details are available in:

- [`docs/extended_readme.md`](docs/extended_readme.md) — longer project description,
- [`docs/workflow_overview.md`](docs/workflow_overview.md) — workflow steps and file flow,
- [`docs/configuration.md`](docs/configuration.md) — configuration options,
