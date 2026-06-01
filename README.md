# Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes

This repository contains a project for the **Modelling of Complex Biological Systems** course.

The project investigates whether inflammatory bowel disease (IBD) is better reflected by changes in the **taxonomic composition** of the gut microbiome, changes in its **functional metabolic potential**, or a combination of both.

The core idea is:

> IBD-related dysbiosis may not only be visible in *which bacteria are present*, but also in *what metabolic functions the microbial community can perform*.

An additional exploratory part of the project compares simple microbial co-abundance networks between healthy and IBD samples.

---

## Background

Inflammatory bowel disease (IBD) is a chronic inflammatory condition of the gastrointestinal tract that includes **Crohn’s disease** and **ulcerative colitis**. Its exact cause is not fully understood, but it is associated with immune dysregulation, environmental factors, and changes in the gut microbiome.

The gut microbiome can be described at two complementary levels:

- **Taxonomic composition** — which microbial taxa are present and in what abundance.
- **Functional composition** — which metabolic pathways are encoded by the microbial community.

Functional profiles are important because different bacterial taxa can sometimes perform similar biological roles. Therefore, disease-related changes may be visible not only at the level of individual species, but also at the level of microbial metabolic potential.

In this project, functional composition is represented by **metabolic pathways**, for example pathways involved in carbohydrate metabolism, amino acid biosynthesis, short-chain fatty acid metabolism, bile acid metabolism, or bacterial cell wall biosynthesis.

---

## Research Question

The main question of this project is:

> Do functional pathway profiles contain comparable or stronger IBD-related signal than taxonomic profiles?

More specifically, the project asks:

1. Do healthy and IBD samples differ in taxonomic diversity?
2. Do healthy and IBD samples differ in functional pathway diversity?
3. Do taxonomic or functional profiles better separate healthy and IBD samples in ordination analysis?
4. Can IBD status be predicted better from taxa, pathways, or both combined?
5. As an exploratory extension, does IBD change the structure of simple microbial co-abundance networks?

---

## Data

The project uses publicly available metagenomic data from:

- `curatedMetagenomicData`

The analysis focuses on stool samples with available disease status.

Samples are divided into two main groups:

- healthy controls,
- IBD patients.

Each sample is represented by two feature matrices:

### 1. Taxonomic profiles

Species-level or genus-level microbial abundance profiles.

Example features:


- *Faecalibacterium prausnitzii*
- *Bacteroides vulgatus*
- *Escherichia coli*
- *Roseburia intestinalis*



### 2. Functional profiles

Metabolic pathway abundance profiles.

Example feature types:

- carbohydrate metabolism pathways
- amino acid biosynthesis pathways
- short-chain fatty acid metabolism pathways
- bile acid metabolism pathways
- cell wall biosynthesis pathways
  
---

## Analysis Overview

The workflow is divided into six main parts:


1. data download
2. preprocessing
3. diversity analysis
4. ordination
5. differential abundance
6. classification
7. exploratory network analysis

## Repository Structure:

```txt
ibd-functional-dysbiosis/
├── README.md
├── Snakefile
├── config/
│   └── config.yaml
├── data/
│   ├── raw/
│   ├── processed/
│   └── metadata/
├── scripts/
│   ├── 01_download_data.R
│   ├── 02_preprocess.R
│   ├── 03_diversity_analysis.py
│   ├── 04_ordination.R
│   ├── 05_differential_abundance.py
│   ├── 06_classification.py
│   └── 07_network_analysis.py
├── notebooks/
│   └── exploratory_analysis.ipynb
├── results/
│   ├── figures/
│   ├── tables/
│   └── models/
├── report/
│   ├── main.tex
│   └── references.bib
├── presentation/
│   └── slides.pdf
└── environment.yml
```

---

## Workflow

The analysis is intended to be reproducible using Snakemake.

To run the full workflow:

`snakemake --cores 4`

Main expected workflow outputs:

```txt
data/processed/taxa_matrix.csv
data/processed/pathway_matrix.csv
data/processed/metadata.csv
results/tables/classification_metrics.csv
results/tables/differential_taxa.csv
results/tables/differential_pathways.csv
results/tables/network_metrics.csv
results/figures/
```

### 1. Data Preprocessing

The first step prepares matched taxonomic, functional, and metadata tables.

Main preprocessing steps:

download selected IBD-related metagenomic data,
select stool samples,
keep samples with known disease status,
define two groups: healthy and IBD,
match sample IDs between metadata, taxa, and pathway matrices,
remove rare taxa and rare pathways,
transform abundance values using log1p or CLR transformation,
save processed matrices for downstream analysis.

Expected outputs:
```
data/processed/taxa_matrix.csv
data/processed/pathway_matrix.csv
data/processed/metadata.csv
```

### 2. Diversity Analysis

This step checks whether IBD is associated with altered microbiome diversity.

At the taxonomic level:

- richness,
- Shannon diversity.

At the functional level:

- pathway richness,
- pathway Shannon diversity.

The main comparison is:

`healthy vs IBD`

Statistical testing:

- Wilcoxon rank-sum test,
- optional FDR correction.

Expected outputs:

```txt
results/figures/diversity_taxa.png
results/figures/diversity_pathways.png
results/tables/diversity_results.csv
```

### 3. Ordination Analysis

Ordination is used to visualize high-dimensional microbiome profiles in two dimensions.

Planned methods:

- PCA for transformed taxonomic profiles,
- PCA for transformed pathway profiles,
- optional PCoA using Bray-Curtis distances.

This analysis asks whether healthy and IBD samples separate better based on:

- taxonomic profiles,
- functional pathway profiles.

If possible, PERMANOVA will be used to test whether disease status explains a significant part of between-sample variation.

Expected outputs:

```txt
results/figures/pca_taxa.png
results/figures/pca_pathways.png
results/tables/permanova_results.csv
```

### 4. Differential Abundance Analysis

This step identifies individual taxa and pathways associated with IBD status.

Planned analysis:

- Wilcoxon rank-sum test for each taxon,
- Wilcoxon rank-sum test for each pathway,
- Benjamini-Hochberg FDR correction,
- ranking of features by adjusted p-value and effect size.

Expected outputs:

```
results/tables/differential_taxa.csv
results/tables/differential_pathways.csv
results/figures/top_taxa_ibd.png
results/figures/top_pathways_ibd.png
```

### 5. Classification Modelling

Classification is used to compare how much IBD-related information is present in different feature sets.

Three models will be compared:

| Model        | Input features             | Purpose                                                     |
| ------------ | -------------------------- | ----------------------------------------------------------- |
| Taxa-only    | taxonomic abundance matrix | tests signal in microbial composition                       |
| Pathway-only | pathway abundance matrix   | tests signal in functional metabolic potential              |
| Combined     | taxa + pathways            | tests whether both levels provide complementary information |

Planned algorithms:

- logistic regression,
- random forest, if time allows.

Evaluation:

- stratified cross-validation,
- ROC-AUC,
- balanced accuracy,
- sensitivity,
- specificity.

Expected outputs:

```
results/tables/classification_metrics.csv
results/figures/roc_auc_comparison.png
results/tables/feature_importance.csv
```

Interpretation:

- if the taxa-only model performs best, species composition contains stronger IBD-related signal;
- if the pathway-only model performs best, functional profiles contain stronger IBD-related signal;
- if the combined model performs best, taxa and pathways provide complementary information;
- if all models perform poorly, IBD-associated microbiome changes may be heterogeneous or confounded by other factors.



### 6. Exploratory Network Analysis

As a compact network-level extension, the project includes a simple co-abundance network analysis.

This part is exploratory and will focus on a limited number of features, for example the top 30 most variable taxa or pathways.

Network construction:

- nodes: selected taxa or pathways,
- edges: strong Spearman correlations,
- separate networks for healthy and IBD samples.

Basic network metrics:

- number of nodes,
- number of edges,
- density,
- average degree,
- clustering coefficient.

Expected outputs:

```txt
results/figures/network_healthy.png
results/figures/network_ibd.png
results/tables/network_metrics.csv
```

Main interpretation:

> The network analysis checks whether IBD is associated with altered organization of microbial co-abundance relationships.

This analysis is not intended to fully reconstruct microbial interactions. It is used as a simple systems-level summary of microbiome organization.

---

## Expected Results

The project is expected to produce a compact comparison of IBD-associated dysbiosis at three levels:

1. Taxonomic level
   1. changes in microbial composition and diversity.
2. Functional level
   1. changes in metabolic pathway composition and diversity.
3. Exploratory network level
   1. changes in simple co-abundance structure.

Possible outcomes:

- functional profiles classify IBD as well as or better than taxonomic profiles,
- taxonomic profiles contain stronger disease-related signal,
- combined taxa + pathway profiles improve classification,
- IBD samples show altered diversity but weak classification performance,
- network structure differs between healthy and IBD samples.


---

### Main limitations:

- only one or a small number of datasets will be used,
- disease status is simplified to healthy vs IBD,
- Crohn’s disease and ulcerative colitis may not be analysed separately,
- longitudinal disease activity is not modelled,
- batch effects may only be inspected, not fully corrected,
- network analysis is exploratory and based on correlation, not direct biological interactions.


### Future Work

- separate analysis of Crohn’s disease and ulcerative colitis,
- longitudinal analysis of disease activity,
- batch-aware modelling across multiple studies,
- MaAsLin2-based association testing,
- MOFA-based integration of taxonomic and functional profiles,
- more robust microbiome network inference using SPIEC-EASI or SparCC,
- taxa-pathway bipartite network analysis.


### Working Hypothesis

> *Functional pathway profiles contain IBD-related information that is comparable to or complementary with taxonomic microbiome profiles.*

The broader interpretation is that IBD-associated dysbiosis may involve not only changes in bacterial composition, but also changes in metabolic potential and microbial system organization.