# Configuration

The main configuration file is:

```text
config/config.yaml
```

The file is checked against:

```text
schemas/config.schema.yaml
```

## Paths

```yaml
paths:
  data_raw: "data/raw"
  data_processed: "data/processed"
  results: "results"
  logs: "results/logs"
  figures: "results/figures"
  tables: "results/tables"
```

These paths define where input, intermediate and output files are written.

## Dataset settings

```yaml
dataset:
  study_name: "HMP_2019_ibdmdb"
  body_site: "stool"
  disease_column: "study_condition"
```

Meaning:

| Key              | Meaning                                      |
| ---------------- | -------------------------------------------- |
| `study_name`     | study selected from `curatedMetagenomicData` |
| `body_site`      | body site kept for analysis                  |
| `disease_column` | metadata column with disease labels          |

More details aobut [curatedMetagenomicData](https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html?utm_source=chatgpt.com)

## Metadata settings

```yaml
metadata:
  sample_id_column: "sample_id"
  disease_column: "study_condition"
  healthy_label: "control"
  ibd_label: "IBD"
```

The workflow maps original labels to two groups:

- `control` becomes `healthy`,
- `IBD` stays `IBD`.

## Preprocessing settings

```yaml
preprocessing:
  min_prevalence: 0.1
  transformation: "log1p"
```

`min_prevalence: 0.1` means that a feature must be present in at least 10% of samples.



## Diversity settings

```yaml
diversity:
  calculate_chao1: false
  metrics:
    - richness
    - shannon
    - simpson
    - invsimpson
    - chao1
```

Chao1 is listed in the config, but it is disabled by default because the data are relative abundance profiles.


## Differential abundance settings

```yaml
differential_abundance:
  p_adjust_method: "BH"
  ntop_features: 20
```

`BH` means Benjamini-Hochberg FDR correction.

## Classification settings

```yaml
classification:
  seed: 123
  test_fraction: 0.25
  positive_class: "IBD"
  ntop_features: 20
```

These values control reproducibility, train/test split size and the positive class for ROC-AUC.

## Plotting settings

```yaml
plotting:
  width: 8
  height: 6
  dpi: 300
```

These values control output figure size and resolution.
