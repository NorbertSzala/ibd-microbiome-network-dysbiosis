"""
Which taxa and fucntional pathways are differentially abundant between healthy and IBD samples?
"""

rule differential_abundance:
    input:
        taxa = f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways = f"{DATA_PROCESSED}/pathway_matrix.csv",
        metadata = f"{DATA_PROCESSED}/metadata.csv"
    
    output:
        taxa_results = f"{TABLES}/differential_taxa.csv",
        pathways_results = f"{TABLES}/differential_pathways.csv",
        top_features = f"{TABLES}/differential_top_features.csv",
        pathways_plot = f"{FIGURES}/top_differential_pathways.png",
        taxa_plot = f"{FIGURES}/top_differential_taxa.png",
        pathways_plot_prevalence = f"{FIGURES}/top_differential_pathways_prevalence.png",
        taxa_plot_prevalence = f"{FIGURES}/top_differential_taxa_prevalence.png"


    params:
        seed=config["classification"]["seed"],
        p_adjust_method=config["differential_abundance"]["p_adjust_method"],
        ntop_features=config["differential_abundance"]["ntop_features"],
        plot_width=config["plotting"]["width"],
        plot_height=config["plotting"]["height"],
        dpi=config["plotting"]["dpi"],
        max_label_words_taxa=config["plotting"]["max_label_words_taxa"],
        max_label_words_pathways=config["plotting"]["max_label_words_pathways"],
            
    conda:
        "../envs/r_analysis.yaml"

    log:
        f"{LOGS}/06_differential_abundance.smk"

    message:
        "Running 06 differential_abundance"

    script:
        "../scripts/06_differential_abundance.R"