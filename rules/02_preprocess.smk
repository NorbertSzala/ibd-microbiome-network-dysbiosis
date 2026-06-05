rule preprocess:

    input:
        taxa=f"{DATA_RAW}/taxa.rds",
        pathways=f"{DATA_RAW}/pathways.rds",
        metadata=f"{DATA_RAW}/metadata.csv"

    output:
        taxa=f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways=f"{DATA_PROCESSED}/pathway_matrix.csv",
        metadata=f"{DATA_PROCESSED}/metadata.csv",
        summary=f"{DATA_PROCESSED}/preprocessing_summary.csv"

    params:
        min_prevalence=config["preprocessing"]["min_prevalence"],
        transformation=config["preprocessing"]["transformation"]

    conda:
        "envs/r.yaml"

    log:
        f"{LOGS}/02_preprocess.log"
        
    script:
        "scripts/02_preprocess.R"