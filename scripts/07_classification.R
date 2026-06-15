#!/usr/bin/env Rscript

# ==============================================================================
# Script: 07_classification.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 15-06-2026
#
# Description:
#   Compares the ability of taxonomic and functional pathway profiles to classify
#   samples as healthy or IBD using random forest models.
#
# Inputs:
#   - data/processed/taxa_matrix.csv
#   - data/processed/pathways_matrix.csv
#   - data/processed/metadata.csv
#
# Outputs:
#   - classification_metrics.csv
#   - classification_feature_importance.csv
#   - classification_auc_comparison.png
#   - classification_feature_importance_taxa.png
#   - classification_feature_importance_pathways.png
# ==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(randomForest)
  library(pROC)
})

source("scripts/functions/00_common_functions.R")
source("scripts/functions/07_classification_functions.R")



# ==============================================================================
# 2. Paths
# ==============================================================================

# When running with Snakemake, paths can be passed automatically.
# For manual runs, define paths here.

if (exists("snakemake")) {
  input_files <- snakemake@input
  output_files <- snakemake@output
  params <- snakemake@params
} else {
  stop("This script should be run through snakemake")
}
set.seed(params$seed)

# Create output directories if needed
walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. Load data
# ==============================================================================
message("Loading input data")
taxa_mat <- load_abundance_matrix(input_files$taxa)
pathways_mat <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)


# ==============================================================================
# 2. Prepare classification datasets
# ==============================================================================
message("Preparing classification datasets")

taxa_df <- prepare_classification_data(
  mat = taxa_mat,
  metadata = metadata
)

pathways_df <- prepare_classification_data(
  mat = pathways_mat,
  metadata  = metadata
)

if (!identical(taxa_df$sample_id, pathways_df$sample_id)) {
  stop("Taxa and pathways datasets do not have the same sample order.")
}
# ==============================================================================
# 3. Split dataset to trainin and test subdatasets
# ==============================================================================
message("Spliting dataset into train and test sets")
taxa_split <- split_train_test(
  df = taxa_df,
  test_fraction = params$test_fraction,
  seed = params$seed
)

pathways_split <- split_train_test(
  df = pathways_df,
  test_fraction = params$test_fraction,
  seed = params$seed
)


# ==============================================================================
# 4. Train models
# ==============================================================================
message("Training random forrest models")
taxa_rf <- train_random_forest(
  df = taxa_split$train,
  seed = params$seed
)

pathways_rf <- train_random_forest(
  df = pathways_split$train,
  seed = params$seed
)

# ==============================================================================
# 5. Evaluate models
# ==============================================================================
message("Evaluating models")

taxa_metrics <- evaluate_random_forest(
  model = taxa_rf,
  df = taxa_split$test,
  feature_set = "taxa",
  positive_class = params$positive_class
)

pathways_metrics <- evaluate_random_forest(
  model = pathways_rf,
  df = pathways_split$test,
  feature_set = "pathways",
  positive_class = params$positive_class
)

metrics_df <- bind_rows(taxa_metrics, pathways_metrics)
# ==============================================================================
# 6. Feature imporance
# ==============================================================================
message("Extracting feature importance")
taxa_results <- run_classification_for_feature_set(
  df = taxa_df,
  feature_set = "taxa",
  models = params$models,
  positive_class = params$positive_class,
  test_fraction = params$test_fraction,
  seed = params$seed
)

pathways_results <- run_classification_for_feature_set(
  df = pathways_df,
  feature_set = "pathways",
  models = params$models,
  positive_class = params$positive_class,
  test_fraction = params$test_fraction,
  seed = params$seed
)

metrics_df <- bind_rows(
  taxa_results$metrics,
  pathways_results$metrics
)

predictions_df <- bind_rows(
  taxa_results$predictions,
  pathways_results$predictions
)

importance_df <- bind_rows(
  taxa_results$importance,
  pathways_results$importance
)

# ==============================================================================
# 7. Plots
# ==============================================================================
message("Creating classification plots")
auc_plot <- plot_auc_comparison(metrics_df)

taxa_importance_plot <- plot_feature_importance(
  importance_df = importance_df,
  selected_feature_set = "taxa",
  ntop = params$ntop_features
)

pathways_importance_plot <- plot_feature_importance(
  importance_df = importance_df,
  selected_feature_set = "pathways",
  ntop = params$ntop_features
)

# ==============================================================================
# 8. Save outputs
# ==============================================================================
message("Saving outputs...")
readr::write_csv(metrics_df, output_files$metrics)
readr::write_csv(predictions_df, output_files$predictions)
readr::write_csv(importance_df, output_files$feature_importance)

auc_plot <- plot_model_auc_comparison(metrics_df)

ggsave(
  filename = output_files$auc_plot,
  plot = auc_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

taxa_rf_plot <- plot_classification_features(
  importance_df = importance_df,
  selected_feature_set = "taxa",
  selected_model = "random_forest",
  ntop = params$ntop_features
)

pathways_rf_plot <- plot_classification_features(
  importance_df = importance_df,
  selected_feature_set = "pathways",
  selected_model = "random_forest",
  ntop = params$ntop_features
)

taxa_en_plot <- plot_classification_features(
  importance_df = importance_df,
  selected_feature_set = "taxa",
  selected_model = "elastic_net",
  ntop = params$ntop_features
)

pathways_en_plot <- plot_classification_features(
  importance_df = importance_df,
  selected_feature_set = "pathways",
  selected_model = "elastic_net",
  ntop = params$ntop_features
)

ggsave(
  filename = output_files$taxa_rf_importance_plot,
  plot = taxa_rf_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$pathways_rf_importance_plot,
  plot = pathways_rf_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$taxa_en_importance_plot,
  plot = taxa_en_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$pathways_en_importance_plot,
  plot = pathways_en_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

# ==============================================================================
# 7. Session info
# ==============================================================================

message("Finished successfully.")
sessionInfo()