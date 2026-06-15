"""
Counts Richness, Shannon and Simpson metrics to define diversity across samples    
"""
rule diversity_analysis:
    input:
        taxa=f"{DATA_PROCESSED}/taxa_matrix_raw_filtered.csv",
        pathways=f"{DATA_PROCESSED}/pathway_matrix_raw_filtered.csv",
        metadata=f"{DATA_PROCESSED}/metadata.csv"
    output:
        diversity_table=f"{TABLES}/diversity_results.csv",
        test_table=f"{TABLES}/diversity_tests.csv",
        taxa_plot=f"{FIGURES}/diversity_taxa.png",
        pathway_plot=f"{FIGURES}/diversity_pathways.png"

    params:
        seed=config["classification"]["seed"],
        calculate_chao1=config["diversity"]["calculate_chao1"],
        metrics=config["diversity"]["metrics"],
        plot_width=config["plotting"]["width"],
        plot_height=config["plotting"]["height"],
        dpi=config["plotting"]["dpi"],
        
    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/04_diversity_analysis.log"

    message:
        "Running 04_diversity_analysis.smk"

    script:
        "../scripts/04_diversity_analysis.R"