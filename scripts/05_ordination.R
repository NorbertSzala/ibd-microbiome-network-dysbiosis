#!/usr/bin/env Rscript

# ==============================================================================
# Script: 05_ordination.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 2026-06-06
#
# Description:
#   Performs ordination analysis for taxonomic and functional pathway profiles. The analysis asks whether healthy and IBD samples differ in global microbiome structure.
#
#   Main analyses:
#   - PCA on transformed abundance matrices. - Does main axes split healthy and IBD?
#   - PCoA on Bray-Curtis distance matrices. - are samples closer to each other in groups than between groups?
#   - PERMANOVA to test global group differences. - Does disease status explain significant part of global difference in microbiome expression profile?
#
# Inputs:
#   - data/processed/taxa_matrix.csv
#   - data/processed/pathway_matrix.csv
#   - data/processed/metadata.csv
#
# Outputs:
#   - results/figures/pca_taxa.png
#   - results/figures/pca_pathways.png
#   - results/figures/pcoa_taxa_bray.png
#   - results/figures/pcoa_pathways_bray.png
#   - results/tables/permanova_results.csv
#   - results/tables/ordination_variance_explained.csv
# ==============================================================================



# ==============================================================================
# 1. Setup
# ==============================================================================
suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(ggplot2)
    library(vegan)
    library(stats)
})

source("scripts/functions/00_common_functions.R")
source("scripts/functions/05_ordination_functions.R")



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
walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)
set.seed(params$seed)

# ==============================================================================
# 3. Load data
# ==============================================================================
message("Loading input data")
taxa_mat <- load_abundance_matrix(input_files$taxa)
pathway_mat <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)


# ==============================================================================
# 4. Match samples
# ==============================================================================
message("Matching samples")
taxa_matched <- match_samples(taxa_mat, metadata)
pathway_matched <- match_samples(pathway_mat, metadata)

taxa_mat <- taxa_matched$mat
taxa_metadata <- taxa_matched$metadata

pathway_mat <- pathway_matched$mat
pathway_metadata <- pathway_matched$metadata


# ==============================================================================
# 5. PCA
# ==============================================================================
message("Running PCA")

pca_taxa <- run_pca(
  mat = taxa_mat,
  metadata = taxa_metadata,
  feature_set = "taxa",
  pca_scale = params$pca_scale
)

pca_pathways <- run_pca(
  mat = pathway_mat,
  metadata = pathway_metadata,
  feature_set = "pathways",
  pca_scale = params$pca_scale
)

pca_taxa_plot <- plot_pca(
  pca_result = pca_taxa,
  title = "PCA of taxonomic profiles"
)

pca_pathway_plot <- plot_pca(
  pca_result = pca_pathways,
  title = "PCA of pathway profiles"
)


# ==============================================================================
# 6. PCoA
# ==============================================================================
message("Running PCoA on Bray-Curtis distances")

pcoa_taxa <- run_pcoa_bray(
  mat = taxa_mat,
  metadata = taxa_metadata,
  feature_set = "taxa"
)

pcoa_pathways <- run_pcoa_bray(
  mat = pathway_mat,
  metadata = pathway_metadata,
  feature_set = "pathways"
)

pcoa_taxa_plot <- plot_pcoa(
  pcoa_result = pcoa_taxa,
  title = "PCoA of taxonomic profiles using Bray-Curtis distance"
)

pcoa_pathway_plot <- plot_pcoa(
  pcoa_result = pcoa_pathways,
  title = "PCoA of pathway profiles using Bray-Curtis distance"
)


# ==============================================================================
# 7. PERMANOVA
# ==============================================================================
message("Running PERMANOVA")

permanova_taxa <- run_permanova(
  mat = taxa_mat,
  metadata = taxa_metadata,
  feature_set = "taxa"
)

permanova_pathways <- run_permanova(
  mat = pathway_mat,
  metadata = pathway_metadata,
  feature_set = "pathways"
)

permanova_results <- bind_rows(
  permanova_taxa,
  permanova_pathways
)


# ==============================================================================
# 8. Save outputs
# ==============================================================================

message("Saving outputs")

ordination_variance <- bind_rows(
  pca_taxa$variance,
  pca_pathways$variance,
  pcoa_taxa$variance,
  pcoa_pathways$variance
)

readr::write_csv(permanova_results, output_files$permanova)
readr::write_csv(ordination_variance, output_files$variance)

ggsave(
  filename = output_files$pca_taxa,
  plot = pca_taxa_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$pca_pathways,
  plot = pca_pathway_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$pcoa_taxa,
  plot = pcoa_taxa_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

ggsave(
  filename = output_files$pcoa_pathways,
  plot = pcoa_pathway_plot,
  width = params$plot_width,
  height = params$plot_height,
  dpi = params$dpi
)

message("Done")
sessionInfo()