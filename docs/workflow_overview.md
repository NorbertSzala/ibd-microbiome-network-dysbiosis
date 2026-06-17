# Workflow overview

The pipeline is controlled by `Snakefile` and rule files in `rules/`.

## Run command

```bash
snakemake --cores 4 --use-conda
```

A dry run can be started with:

```bash
snakemake -n --use-conda
```

## Workflow steps

| Step | Rule                     | Script                                | Main task                                                                 |
| ---- | ------------------------ | ------------------------------------- | ------------------------------------------------------------------------- |
| 01   | `download_data`          | `scripts/01_download_data.R`          | Download taxa, pathway and metadata objects from `curatedMetagenomicData` |
| 02   | `preprocess`             | `scripts/02_preprocess.R`             | Match samples, filter features, transform matrices                        |
| 03   | `abundance_distribution` | `scripts/03_abundance_distribution.R` | Check abundance distributions after preprocessing                         |
| 04   | `diversity_analysis`     | `scripts/04_diversity_analysis.R`     | Calculate alpha-diversity and group tests                                 |
| 05   | `ordination`             | `scripts/05_ordination.R`             | Run PCA, PCoA, PERMANOVA and PERMDISP                                     |
| 06   | `differential_abundance` | `scripts/06_differential_abundance.R` | Test taxa and pathways one by one                                         |
| 07   | `classification`         | `scripts/07_classification.R`         | Train and evaluate classifiers                                            |

## File flow

```text
curatedMetagenomicData
        |
        v
01_download_data
        |
        v
data/raw/taxa.rds
data/raw/pathways.rds
data/raw/metadata.csv
        |
        v
02_preprocess
        |
        v
data/processed/taxa_matrix.csv
data/processed/pathway_matrix.csv
data/processed/metadata.csv
        |
        v
03-07 analysis rules
        |
        v
results/tables/
results/figures/
```

## Main final outputs

The main rule creates analysis tables and figures in:

```text
results/tables/
results/figures/
```

The workflow also writes logs to:

```text
results/logs/
```

## Notes

- Rules 01 and 02 prepare the input matrices.
- Rules 03-07 perform the main analysis.
- Conda environments are defined in `envs/`.
- Configuration is validated with `schemas/config.schema.yaml`.
