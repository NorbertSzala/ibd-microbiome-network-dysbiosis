# snakemake --cores 4 --use-conda --conda-frontend mamba

# ------------------------------------------------------------------------------
# --- Imports ------------------------------------------------------------------
# ------------------------------------------------------------------------------

configfile: "config/config.yaml"

import pandas as pd
from glob import glob
from snakemake.utils import validate
from pathlib import Path

validate(configfile, "schemas/config.schema.yaml")


# ------------------------------------------------------------------------------
# --- Path variables -----------------------------------------------------------
# ------------------------------------------------------------------------------

DATA_RAW = CONFIG['paths']["data_raw"]
DATA_METADATA = CONFIG['paths']["data_metadata"]
DATA_PROCESSED = CONFIG['paths']["data_processed"]

RESULTS = CONFIG['paths']["results"]
LOGS = CONFIG['paths']["results/logs"]

FIGURES = CONFIG['paths']["results/figures"]
MODELS = CONFIG['paths']["results/models"]
TABLES = CONFIG['paths']["results/tables"]


# ------------------------------------------------------------------------------
# --- Helper functions ---------------------------------------------------------
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# --- Main Rule ----------------------------------------------------------------
# ------------------------------------------------------------------------------

rule all:
    input:
        expand(
            
        )

# ------------------------------------------------------------------------------
# --- Include Rules ------------------------------------------------------------
# ------------------------------------------------------------------------------

include "rules/01_download_data.smk"
include "rules/02_preprocess.smk"
include "rules/03_diversity_analysis.smk"
include "rules/04_ordination.smk"
include "rules/05_differential_abundance.smk"
include "rules/06_classification.smk"
include "rules/07_network_analysis.smk"