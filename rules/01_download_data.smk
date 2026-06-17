
"""
Description:
    Downloads taxonomic relative abundance profiles and functional pathway
    abundance profiles from curatedMetagenomicData. The script saves raw
    SummarizedExperiment objects and sample metadata for downstream preprocessing.

Inputs:
    No file inputs.

Outputs:
    - data/raw/taxa.rds
    - data/raw/pathways.rds
    - data/raw/metadata.csv
    - data/raw/sample_summary.csv
"""

rule download_data:
    output:
        taxa = f"{DATA_RAW}/taxa.rds",
        pathways = f"{DATA_RAW}/pathways.rds",
        metadata = f"{DATA_RAW}/metadata.csv",
        summary = f"{DATA_RAW}/sample_summary.csv"
        
    params:
        study_name = config['dataset']['study_name'],
        body_site = config['dataset']['body_site'],
        disease_column = config['dataset']['disease_column']


    conda:
        "../envs/r_download_data.yaml"

    log:
        f"{LOGS}/01_download_data.log"

    script:
        "../scripts/01_download_data.R"
