"""
Create QC plots showing the distribution of transformed abundance values
after preprocessing.
"""
    
rule abundance_distribution:
    input:
        taxa = f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways = f"{DATA_PROCESSED}/pathway_matrix.csv",

    output:
        summary_table = f"{TABLES}/abundance_distribution_summary.csv",
        taxa_ecdf_plot = f"{FIGURES}/abundance_distribution_ecdf_taxa.png",
        pathways_ecdf_plot = f"{FIGURES}/abundance_distribution_ecdf_pathways.png",
        taxa_hist_plot = f"{FIGURES}/abundance_distribution_hist_taxa.png",
        pathways_hist_plot = f"{FIGURES}/abundance_distribution_hist_pathways.png", 

    params:
        plot_width=config["plotting"]["width"],
        plot_height=config["plotting"]["height"],
        dpi=config["plotting"]["dpi"],

    conda:
        "../envs/r_analysis.yaml"

    log:
        f"{LOGS}/03_abundance_distribution.log"

    message:
        "Runnign 03_abundance_distribution rule"

    script:
        "../scripts/03_abundance_distribution.R"

