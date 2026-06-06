"""
Simple rule making plot about abundance distribution samples after processing    
"""
    
rule abundance_distribution:
    input:
        taxa = f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways = f"{DATA_PROCESSED}/pathway_matrix.csv",
        metadata = f"{DATA_PROCESSED}/metadata.csv"

    output:
        summary_table = f"{TABLES}/abundance_distribution_summary.csv",
        taxa_density_plot = f"{FIGURES}/abundance_distribution_dens_taxa.png",
        pathways_density_plot = f"{FIGURES}/abundance_distribution_dens_pathways.png",
        taxa_hist_plot = f"{FIGURES}/abundance_distribution_hist_taxa.png",
        pathways_hist_plot = f"{FIGURES}/abundance_distribution_hist_pathways.png", 

    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/03_abundance_distribution.log"

    message:
        "Runnign 03_abundance_distribution rule"

    script:
        "../scripts/03_abundance_distribution.R"

