rule classification:
    input:
        taxa = f"{DATA_PROCESSED}/taxa_matrix.csv",
        pathway_matrix = f"{DATA_PROCESSED}/pathway_matrix.csv",
        metadata= f"{DATA_PROCESSED}/metadata.csv"
        
    output:
        metrics=f"{TABLES}/classification_metrics.csv",
        predictions=f"{TABLES}/classification_predictions.csv",
        feature_importance=f"{TABLES}/classification_feature_importance.csv",
        auc_plot=f"{FIGURES}/classification_auc_comparison.png",
        taxa_rf_importance_plot=f"{FIGURES}/classification_rf_importance_taxa.png",
        pathways_rf_importance_plot=f"{FIGURES}/classification_rf_importance_pathways.png",
        taxa_en_importance_plot=f"{FIGURES}/classification_elastic_net_coefficients_taxa.png",
        pathways_en_importance_plot=f"{FIGURES}/classification_elastic_net_coefficients_pathways.png"

    params:
        seed=config["classification"]["seed"],
        test_fraction=config["classification"]["test_fraction"],
        positive_class=config["classification"]["positive_class"],
        ntop_features=config["classification"]["ntop_features"],
        plot_width=config["plotting"]["width"],
        plot_height=config["plotting"]["height"],
        dpi=config["plotting"]["dpi"],

        

    conda:
        "../envs/r.yaml"

    log:
        f"{LOGS}/07_classification.smk"

    message:
        "Running classification rule."

    script:
        "../scripts/07_classification.R"
