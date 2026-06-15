#!/usr/bin/env Rscript

# ==============================================================================
# Script: 04_diversity_analysis.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 2026-06-06
#
# Description:
#   This script calculates alpha diversity for taxonomic and functional pathway abundance profiles. Alpha diversity describes diversity inside a single sample.
#   The analysis asks whether IBD samples have lower or different taxa/pathway  diversity compared with healthy controls.
#
#   The input matrices contain relative abundance values, not integer raw counts. Because of that, Chao1 is not calculated here. Chao1 requires raw count data, because it depends on singleton and doubleton counts. For normalized abundance profiles this would not be methodologically correct.
#
#   Metrics calculated in this script:
#   - Richness:
#       Number of detected taxa/pathways in a sample. This counts all features
#       with abundance greater than zero.
#
#   - Thresholded richness:
#       Number of taxa/pathways with abundance greater than 1e-5. This is useful
#       for relative abundance data because very small non-zero values may be
#       unstable or biologically weak signals.
#
#   - Shannon index:
#       Measures diversity by combining richness and evenness. Higher Shannon
#       values mean that it is harder to predict the identity of a randomly chosen
#       feature from the sample.
#
#   - Simpson index:
#       Measures dominance structure. It is calculated as 1 - sum(p_i^2).
#       Higher values indicate higher diversity and lower dominance by only a few
#       abundant features.
#
#   - Inverse Simpson index:
#       Gives the effective number of dominant features. It is easier to interpret
#       than the raw Simpson value.
#
#   - Hill numbers:
#       Hill q0 is richness, Hill q1 is exp(Shannon), and Hill q2 is inverse
#       Simpson. They describe the effective number of features at different
#       sensitivity levels to rare and dominant features.
#
#   - Pielou evenness:
#       Measures how evenly abundance is distributed across detected features.
#       Values closer to 1 mean a more even community or pathway profile.
#


# Worth to read:
# https://biostatsquid.com/alpha-diversity-metrics/
# https://doi.org/10.1038/s41598-024-77864-y

# Inputs:
#   - data/processed/taxa_matrix.csv
#   - data/processed/pathways_matrix.csv
#   - data/processed/metadata.csv
#
# Outputs:
#   - results/tables/diversity_results.csv
#   - results/tables/diversity_tests.csv
#   - results/figures/diversity_taxa.png
#   - results/figures/diversity_pathways.png
#
# Notes:
#   Chao1 is not calculated because the available matrices contain relative abundance values. If raw integer count matrices become available, Chao1 can be calculated in a separate count-based diversity analysis.==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(vegan)
})
source("scripts/functions/00_common_functions.R")
source("scripts/functions/04_diversity_analysis_functions.R")


# ==============================================================================
# 2. Paths
# ==============================================================================

# When running with Snakemake, paths can be passed automatically.
# For manual runs, define paths here.

if (exists("snakemake")) {
  input_files <- snakemake@input
  output_files <- snakemake@output
  params<- snakemake@params
} else {
    stop("This script should be run through Snakemake")
}

walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)

set.seed(params$seed)


# ==============================================================================
# 3. Load data
# ==============================================================================

message("Loading input data...")
taxa_mat <- load_abundance_matrix(input_files$taxa)
pathways_mat <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)


# ==============================================================================
# 4. Calculate diversity metrics
# ==============================================================================
message("Calculating diversity metrics")

taxa_diversity <- calculate_diversity_metrics_vegan(
  mat = taxa_mat,
  feature_set = "taxa",
  calculate_chao1 = FALSE
)

pathways_diversity <- calculate_diversity_metrics_vegan(
  mat = pathways_mat,
  feature_set = "pathways",
  calculate_chao1 = FALSE
)

diversity_results <- bind_rows(
  taxa_diversity,
  pathways_diversity
) %>%
  join_metadata(metadata)



# ==============================================================================
# 5. Statistical tests
# ==============================================================================
message("Running Wilcoxon tests")
diversity_tests <- run_wilcoxon_tests(diversity_results)


# ==============================================================================
# 6. Create plots
# ==============================================================================
message("Creating diversity plots")

taxa_plot <- create_diversity_plot(
  diversity_df = diversity_results,
  selected_feature_set = "taxa",
  title = "Taxonomic alpha diversity in healthy and IBD samples"
)

pathways_plot <- create_diversity_plot(
  diversity_df = diversity_results,
  selected_feature_set = "pathways",
  title = "Functional pathway alpha diversity in healthy and IBD samples"
)

# ==============================================================================
#7. Save outputs
# ==============================================================================
message("Saving outputs...")

readr::write_csv(diversity_results, output_files$diversity_table)
readr::write_csv(diversity_tests, output_files$test_table)

ggsave(
  filename = output_files$taxa_plot,
  plot = taxa_plot,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)

ggsave(
  filename = output_files$pathway_plot,
  plot = pathways_plot,
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)


message("Finished successfully.")
sessionInfo()