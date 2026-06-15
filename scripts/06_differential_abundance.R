#!/usr/bin/env Rscript

# ==============================================================================
# Script: 06_differential_abundance.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 2026-06-13
#
# Description:
#   Which taxa and fucntional pathways are differentially abundant between healthy and IBD samples?
#
# Inputs:
#   data/processed/taxa_matrix.csv
#   data/processed/pathways_matrix.csv
#   data/processed/metadata.csv
#
# Outputs:
#   differential_taxa.csv
#   differential_pathways.csv
#   differential_top_features.csv
#   top_differential_pathways.png
#   top_differential_taxa.png
#
# Notes:
#   This script is intended to be run directly or through Snakemake.
# ==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================
source("scripts/functions/00_common_functions.R")
source("scripts/functions/06_differential_abundance_functions.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
})


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
    stop("This script should be run through Snakemake")
}

# Create output directories if needed
walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)
set.seed(params$seed)

# ==============================================================================
# 3. Load data
# ==============================================================================
message("Loading input data")
taxa_mat <- load_abundance_matrix(input_files$taxa)
pathways_mat <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)

# ==============================================================================
# 4. Match samples
# ==============================================================================
message("Match samples")
taxa_matched <- match_samples(taxa_mat, metadata)
pathways_matched <- match_samples(pathways_mat, metadata)

taxa_mat <- taxa_matched$mat
taxa_metadata <- taxa_matched$metadata

pathways_mat <- pathways_matched$mat
pathways_metadata <- pathways_matched$metadata

# ==============================================================================
# 5. Transform data to long format
# ==============================================================================
message("Convert matrivs to long df format")
taxa_long_df <- matrix_to_long_df(mat = taxa_mat, feature_set = "taxa") %>% 
  left_join(taxa_metadata %>% 
    select(sample_id, disease_status), 
    by = "sample_id")

pathways_long_df <- matrix_to_long_df(mat = pathways_mat, feature_set = "pathways") %>% 
  left_join(pathways_metadata %>% 
    select(sample_id, disease_status),
    by = "sample_id")

if (any(is.na(taxa_long_df$disease_status))) {
  stop("Missing disease_status values in taxa_long_df.")
}

if (any(is.na(pathways_long_df$disease_status))) {
  stop("Missing disease_status values in pathways_long_df.")
}


# ==============================================================================
# 6. Differential abundance analysis
# ==============================================================================
message("Running differential abundance analysis")
taxa_results <- run_differential_abundance(taxa_long_df, p_adjust_method = params$p_adjust_method)
pathways_results <- run_differential_abundance(pathways_long_df, p_adjust_method = params$p_adjust_method)

all_results <- bind_rows(taxa_results, pathways_results)

top_features <- get_top_features(results_df = all_results, ntop = params$ntop_features)



# ==============================================================================
# 7. Make plots
# ==============================================================================
message("Making plots")
taxa_plot <- plot_top_features(
  results_df = taxa_results,
  selected_feature_set = "taxa",
  ntop = params$ntop_features,
  title = paste("Top ", params$ntop_features, " differentially abundand taxa")
)

pathways_plot <- plot_top_features(
  results_df = pathways_results,
  selected_feature_set = "pathways",
  ntop = params$ntop_features,
  title = paste("Top ", params$ntop_features, " differentially abundand pathways")
)


taxa_plot_prevalence <- plot_top_prevalence_features(
  results_df = taxa_results,
  selected_feature_set = "taxa",
  ntop = params$ntop_features,
  max_label_words_pathways = params$max_label_words_pathways,
  max_label_words_taxa = params$max_label_words_taxa,
  title = paste("Top", params$ntop_features, "taxa by prevalence difference")
)

pathways_plot_prevalence <- plot_top_prevalence_features(
  results_df = pathways_results,
  selected_feature_set = "pathways",
  ntop = params$ntop_features,
  max_label_words_pathways = params$max_label_words_pathways,
  max_label_words_taxa = params$max_label_words_taxa,
  title = paste("Top", params$ntop_features, "pathways by prevalence difference")
)


# ==============================================================================
# 8. Save outputs
# ==============================================================================
message("Saving outputs...")
readr::write_csv(taxa_results, output_files$taxa_results)
readr::write_csv(pathways_results, output_files$pathways_results)
readr::write_csv(top_features, output_files$top_features)


ggsave(
  filename = output_files$taxa_plot,
  plot = taxa_plot,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)

ggsave(
  filename = output_files$pathways_plot,
  plot = pathways_plot,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)

ggsave(
  filename = output_files$taxa_plot_prevalence,
  plot = taxa_plot_prevalence,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)

ggsave(
  filename = output_files$pathways_plot_prevalence,
  plot = pathways_plot_prevalence,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)

# ==============================================================================
# 9. Session info
# ==============================================================================

message("Finished successfully.")
sessionInfo()