rule preprocess:
    input:
        taxa=f"{DATA_RAW}/taxa.rds",
        pathways=f"{DATA_RAW}/pathways.rds",
        metadata=f"{DATA_RAW}/metadata.csv"

    output:
        taxa=f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways=f"{DATA_PROCESSED}/pathway_matrix.csv",
        taxa_raw_filtered=f"{DATA_PROCESSED}/taxa_matrix_raw_filtered.csv",
        pathways_raw_filtered=f"{DATA_PROCESSED}/pathway_matrix_raw_filtered.csv",
        metadata=f"{DATA_PROCESSED}/metadata.csv",
        summary=f"{DATA_PROCESSED}/preprocessing_summary.csv"

    params:
        min_prevalence=config["preprocessing"]["min_prevalence"],
        transformation=config["preprocessing"]["transformation"],
        sample_id_column=config["metadata"]["sample_id_column"],
        disease_column=config["metadata"]["disease_column"],
        healthy_label=config["metadata"]["healthy_label"],
        ibd_label=config["metadata"]["ibd_label"]

    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/02_preprocess.log"
        
    script:
        "../scripts/02_preprocess.R"