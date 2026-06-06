"""
Description:
    Performs ordination analysis for taxonomic and pathway profiles.
    The script generates PCA and PCoA plots and runs PERMANOVA to test
    whether healthy and IBD samples differ in global microbiome structure.
"""
rule ordination:
    input:
        taxa=f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathways=f"{DATA_PROCESSED}/pathway_matrix.csv",
        metadata=f"{DATA_PROCESSED}/metadata.csv"

    output:
        pca_taxa=f"{FIGURES}/pca_taxa.png",
        pca_pathways=f"{FIGURES}/pca_pathways.png",
        pcoa_taxa=f"{FIGURES}/pcoa_taxa_bray.png",
        pcoa_pathways=f"{FIGURES}/pcoa_pathways_bray.png",
        permanova=f"{TABLES}/permanova_results.csv",
        variance=f"{TABLES}/ordination_variance_explained.csv"

    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/05_ordination.log"
    
    message:
        "Running 05 ordination"

    script:
        "../scripts/05_ordination.R"