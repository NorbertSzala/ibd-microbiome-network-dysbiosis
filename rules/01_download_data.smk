rule rulename:
    
    """
    # Description:
    #   Downloads taxonomic relative abundance profiles and functional pathway
    #   abundance profiles from curatedMetagenomicData. The script saves raw
    #   SummarizedExperiment objects and sample metadata for downstream preprocessing.
    #
    # Inputs:
    #   No file inputs.
    
    """

    output:
        taxa = f"{DATA_RAW}/taxa.rds"
        pathways = f"{DATA_RAW}/pathways.rds"
        metadata = f"{DATA_RAW}/metadata.csv"
        summary = f"{DATA_RAW}/sample_summary.csv"
        
    params:
        study_name = config['dataset']['study_name']
        body_site = config['dataset']['body_site']

    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/01_download_data.log"

    script:
        """
        .../scripts/01_download_data.R
        """
