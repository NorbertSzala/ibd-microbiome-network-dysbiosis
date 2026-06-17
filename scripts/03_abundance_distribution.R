#!/usr/bin/env Rscript

# ==============================================================================
# Script: 03_abundance_distribution.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 06.04.2026

# Description:
#   Creates basic abundance distribution plots for processed taxonomic and functional pathway matrices. This is a simple quality-control and exploratory visualization step after preprocessing.
#
# Inputs:
#   - data/processed/t  axa_matrix.csv
#   - data/processed/pathways_matrix.csv
#   - data/processed/metadata.csv
#
# Outputs:
#   - results/figures/abundance_distribution_taxa.png
#   - results/figures/abundance_distribution_pathways.png
#   - results/tables/abundance_distribution_summary.csv
# ==============================================================================

# ==============================================================================
# 1. Setup
# ==============================================================================
suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
})

source("scripts/functions/00_common_functions.R")
source("scripts/functions/03_abundance_distribution_functions.R")
# ==============================================================================
# 2. Paths
# ==============================================================================

if (exists("snakemake")) {
    input_files <- snakemake@input
    output_files <- snakemake@output
    params <- snakemake@params
} else {
    stop("This script should be run through Snakemake")
}

walk(dirname(unlist(output_files)),
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
)


# ==============================================================================
# 3. Load data
# ==============================================================================
message("Loading processed matrices")
taxa <- load_abundance_table(input_files$taxa)
pathways <- load_abundance_table(input_files$pathways)

# ==============================================================================
#4. Convert to long format
# ==============================================================================
message("Converting processed data to long format")

taxa_long <- matrix_to_long(taxa, feature_set = 'taxa')
pathways_long <- matrix_to_long(pathways, feature_set = 'pathways')

# ==============================================================================
# 5. Plot abundance distribution
# ==============================================================================
message("Creating plots")
taxa_ecdf_plot <- create_abundance_ecdf_plot(
    taxa_long,
    title = "Cumulative distribution of transformed taxonomic abundances"
)
pathways_ecdf_plot <- create_abundance_ecdf_plot(
    pathways_long,
    title = "Cumulative distribution of transformed pathway abundances"
)


taxa_hist_plot<- create_abundance_hist_plot(
    taxa_long, 
    title = "Distribution of transformed taxonomic abundances"
)

pathways_hist_plot<- create_abundance_hist_plot(
    pathways_long, 
    title = "Distribution of transformed pathway abundances"
)


ggsave(
    filename = output_files$taxa_ecdf_plot, 
    plot = taxa_ecdf_plot, 
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)
ggsave(
    filename = output_files$pathways_ecdf_plot, 
    plot = pathways_ecdf_plot, 
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)


ggsave(
    filename = output_files$taxa_hist_plot, 
    plot = taxa_hist_plot, 
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)
ggsave(
    filename = output_files$pathways_hist_plot, 
    plot = pathways_hist_plot, 
    width = params$plot_width,
    height = params$plot_height,
    dpi = params$dpi
)
# ==============================================================================
# 6. Save summary
# ==============================================================================

message("Saving abundance summary")

summary_df <- bind_rows(
    summarize_abundance(taxa_long),
    summarize_abundance(pathways_long)
)

readr::write_csv(summary_df, output_files$summary_table)
message("Done")
sessionInfo()