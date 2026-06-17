# snakemake --cores 4 --use-conda --conda-frontend mamba

# ------------------------------------------------------------------------------
# --- Imports ------------------------------------------------------------------
# ------------------------------------------------------------------------------

configfile: "config/config.yaml"

import pandas as pd
from glob import glob
from snakemake.utils import validate
from pathlib import Path

validate(config, "schemas/config.schema.yaml")


# ------------------------------------------------------------------------------
# --- Path variables -----------------------------------------------------------
# ------------------------------------------------------------------------------

DATA_RAW = config['paths']["data_raw"]
DATA_METADATA = config['paths']["data_metadata"]
DATA_PROCESSED = config['paths']["data_processed"]

RESULTS = config['paths']["results"]
LOGS = config['paths']["logs"]

FIGURES = config['paths']["figures"]
MODELS = config['paths']["models"]
TABLES = config['paths']["tables"]



# ------------------------------------------------------------------------------
# --- FInal outputs ------------------------------------------------------------
# ------------------------------------------------------------------------------

FINAL_OUTPUTS = [
        # 02_preprocess
        # f"{DATA_PROCESSED}/taxa_matrix.csv",
        # f"{DATA_PROCESSED}/pathway_matrix.csv",
        # f"{DATA_PROCESSED}/metadata.csv",
        # f'{DATA_PROCESSED}/taxa_matrix_raw_filtered.csv',
        # f"{DATA_PROCESSED}/pathway_matrix_raw_filtered.csv",
        # f"{DATA_PROCESSED}/preprocessing_summary.csv",

        # 03_abundance_distribution
        f"{FIGURES}/abundance_distribution_dens_taxa.png",
        f"{FIGURES}/abundance_distribution_dens_pathways.png",
        f"{FIGURES}/abundance_distribution_hist_taxa.png",
        f"{FIGURES}/abundance_distribution_hist_pathways.png", 
        f"{TABLES}/abundance_distribution_summary.csv",
        
        #04_diversity_analysis
        f"{TABLES}/diversity_results.csv",
        f"{TABLES}/diversity_tests.csv",
        f"{FIGURES}/diversity_taxa.png",
        f"{FIGURES}/diversity_pathways.png",

        # 05_ordination
        f"{FIGURES}/pca_taxa.png",
        f"{FIGURES}/pca_pathways.png",
        f"{FIGURES}/pcoa_taxa_bray.png",
        f"{FIGURES}/pcoa_pathways_bray.png",
        f"{TABLES}/permanova_results.csv",
        f"{TABLES}/ordination_variance_explained.csv",


        # 06_differential_abundance
        f"{TABLES}/differential_taxa.csv",
        f"{TABLES}/differential_pathways.csv",
        f"{TABLES}/differential_top_features.csv",
        f"{FIGURES}/top_differential_pathways.png",
        f"{FIGURES}/top_differential_taxa.png",

        # 07_classifications
        f"{TABLES}/classification_metrics.csv",
        f"{TABLES}/classification_predictions.csv",
        f"{TABLES}/classification_feature_importance.csv",
        f"{FIGURES}/classification_auc_comparison.png",
        f"{FIGURES}/classification_rf_importance_taxa.png",
        f"{FIGURES}/classification_rf_importance_pathways.png",
        f"{FIGURES}/classification_elastic_net_coefficients_taxa.png",
        f"{FIGURES}/classification_elastic_net_coefficients_pathways.png"
]

# ------------------------------------------------------------------------------
# --- Main rule --- ------------------------------------------------------------
# ------------------------------------------------------------------------------
rule all:
    input:
        FINAL_OUTPUTS


# ------------------------------------------------------------------------------
# --- Include Rules ------------------------------------------------------------
# ------------------------------------------------------------------------------
# include: "rules/01_download_data.smk"
# include: "rules/02_preprocess.smk"
include: "rules/03_abundance_distribution.smk"
include: "rules/04_diversity_analysis.smk"
include: "rules/05_ordination.smk"
include: "rules/06_differential_abundance.smk"
include: "rules/07_classification.smk"
