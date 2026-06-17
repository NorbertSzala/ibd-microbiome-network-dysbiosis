# Extended README

## Project title

**Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes**

## Aim

The project studies gut microbiome dysbiosis in inflammatory bowel disease (IBD).

The main aim is to compare two types of microbiome information:

- **taxonomic information**: which microbial taxa are present,
- **functional information**: which metabolic pathways are present.

The key question is:

> Do pathway profiles contain IBD-related signal that is similar to, stronger than, or weaker than taxonomic profiles?

## Biological background

IBD is a chronic inflammatory disease of the gastrointestinal tract. It includes Crohn's disease and ulcerative colitis.

IBD is linked with changes in the gut microbiome. These changes can be visible in microbial composition, but also in microbial functions. This is important because different microbes can sometimes perform similar metabolic roles. For this reason, disease signal may be partly taxonomic and partly functional.

In this project, functional profiles are represented by metabolic pathway abundance profiles.

## Input data

The workflow uses public metagenomic profiles from `curatedMetagenomicData`.

The default study is:

```yaml
study_name: "HMP_2019_ibdmdb"
body_site: "stool"
disease_column: "study_condition"
```

The workflow uses two abundance matrices:

1. taxonomic abundance matrix,
2. pathway abundance matrix.

The metadata table is used to select stool samples and define two groups:

- `healthy`,
- `IBD`.

## Analysis design

The same samples are used for taxonomic and pathway analysis when possible. This makes the comparison more direct.

The workflow performs:

1. data download,
2. preprocessing,
3. abundance distribution checks,
4. alpha-diversity analysis,
5. ordination analysis,
6. differential abundance analysis,
7. classification analysis.

## Preprocessing

Preprocessing prepares clean matrices for analysis.

Main steps:

- load taxonomic and pathway profiles,
- load metadata,
- keep stool samples,
- keep samples with known disease status,
- map original labels to `healthy` and `IBD`,
- match samples between taxa, pathways, and metadata,
- remove rare features using a prevalence threshold,
- remove special and taxon-stratified HUMAnN pathway features,
- transform abundance values with `log1p`, unless `none` is selected,
- save processed matrices and metadata.

## Diversity analysis

Diversity analysis checks whether IBD samples have different within-sample diversity.

Calculated metrics include:

- richness,
- Shannon diversity,
- Simpson diversity,
- inverse Simpson diversity.

Chao1 is disabled by default because the available matrices contain relative abundance values, not raw integer counts.

The workflow compares `healthy` and `IBD` samples using Wilcoxon rank-sum tests. P-values are adjusted with the Benjamini-Hochberg method.

## Ordination analysis

Ordination is used to visualise high-dimensional microbiome profiles in two dimensions.

The workflow runs:

- PCA on transformed abundance matrices,
- PCoA on Bray-Curtis distance matrices,
- PERMANOVA to test group-level differences,
- PERMDISP to check whether group dispersion affects PERMANOVA interpretation.

PCA and PCoA are produced separately for taxa and pathways.

## Differential abundance analysis

Differential abundance analysis tests each feature separately.

For each taxon or pathway, the workflow calculates:

- mean abundance in healthy samples,
- mean abundance in IBD samples,
- median abundance in both groups,
- prevalence in both groups,
- Wilcoxon p-value for abundance difference,
- Fisher exact test p-value for prevalence difference,
- adjusted p-values.

Top features are selected using adjusted p-values and effect size.

## Classification analysis

Classification analysis tests how well each feature set predicts disease status.

The workflow compares:

- taxa-only model,
- pathway-only model.

Two model types are used:

- random forest,
- elastic net logistic regression.

The workflow reports metrics such as ROC-AUC, accuracy, balanced accuracy, sensitivity, and specificity.

Feature importance is also saved for both model types.

## Main interpretation

The project should be interpreted as a comparison of signal strength.

Possible outcomes:

- taxa classify IBD better than pathways,
- pathways classify IBD better than taxa,
- both feature sets perform similarly,
- diversity or ordination shows group differences but classification remains weak.

These outcomes help describe whether IBD signal is stronger at the taxonomic or functional level.

