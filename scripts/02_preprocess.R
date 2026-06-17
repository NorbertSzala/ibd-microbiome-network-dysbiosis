#!/usr/bin/env Rscript

# ==============================================================================
# Script: 02_preprocess.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 05-06-2026
#
# Description:
#   Preprocesses raw taxonomic and functional metagenomic profiles downloaded
#   from curatedMetagenomicData.
#
#   Main tasks:
#   - load raw SummarizedExperiment objects(rds format),
#   - extract abundance matrices,
#   - match samples between taxa, pathways, and metadata,
#   - define simplified disease status groups,
#   - keep healthy and IBD samples,
#   - remove rare features,
#   - transform abundance values,
#   - save processed matrices for downstream analysis.
#
# Inputs:
#   - data/raw/taxa.rds
#   - data/raw/pathways.rds
#   - data/raw/metadata.csv
#
# Outputs:
#   - data/processed/taxa_matrix.csv
#   - data/processed/pathway_matrix.csv
#   - data/processed/metadata.csv
#   - data/processed/preprocessing_summary.csv
#
# Notes:
#   This script is intended to be run through Snakemake.
# ==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
  library(readr)
  library(tibble)
  library(stringr)
})

source("scripts/functions/02_preprocess_functions.R")

# ==============================================================================
# 2. Paths and parameters
# ==============================================================================

if (exists("snakemake")) {
  input_files <- snakemake@input
  output_files <- snakemake@output
  params <- snakemake@params
} else {
  stop("This script must be run through Snakemake because input paths, output paths, and parameters are passed by the workflow.")
}
min_prevalence <- as.numeric(params$min_prevalence)
transformation <- as.character(params$transformation)
pseudocount <- as.numeric(params$pseudocount)

sample_id_column <- as.character(params$sample_id_column)
disease_column <- as.character(params$disease_column)
healthy_label <- as.character(params$healthy_label)
ibd_label <- as.character(params$ibd_label)

message("Minimum feature prevalence: ", min_prevalence)
message("Transformation method: ", transformation)
message("Sample ID column: ", sample_id_column)
message("Disease column: ", disease_column)
message("Healthy label: ", healthy_label)
message("IBD label: ", ibd_label)

walk(
  dirname(unlist(output_files)),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Load data
# ==============================================================================
message("Loading raw objects...")

taxa_se <- readRDS(input_files$taxa)
pathway_se <- readRDS(input_files$pathways)
metadata <- readr::read_csv(input_files$metadata, show_col_types = FALSE)

if (!"sample_id" %in% colnames(metadata)) {
  stop("Metadata must contain a `sample_id` column.")
}


# ==============================================================================
# 4. Extract abundance matrices
# ==============================================================================

message("Extracting abundance matrices...")

taxa_mat <- extract_abundance_matrix(taxa_se, assay_name = "relative_abundance")
pathway_mat <- extract_abundance_matrix(pathway_se, assay_name = "pathway_abundance")

message("Raw taxa matrix dimensions: ", paste(dim(taxa_mat), collapse = " x "))
message("Raw pathway matrix dimensions: ", paste(dim(pathway_mat), collapse = " x "))


# ==============================================================================
# 5. Metadata filtering and sample matching
# ==============================================================================

message("Filtering metadata to healthy and IBD samples...")

# Creates standarized labels
metadata <- filter_metadata_to_groups(
  metadata = metadata,
  sample_id_column = sample_id_column,
  disease_column = disease_column,
  healthy_label = healthy_label,
  ibd_label = ibd_label
)

# Extract only those rows common in all datas (kinda innerjoin)
matched <- match_taxa_pathway_samples(
  taxa_mat = taxa_mat,
  pathway_mat = pathway_mat,
  metadata = metadata
)

taxa_mat <- matched$taxa
pathway_mat <- matched$pathways
metadata <- matched$metadata

message("Matched taxa matrix dimensions: ", paste(dim(taxa_mat), collapse = " x "))
message("Matched pathway matrix dimensions: ", paste(dim(pathway_mat), collapse = " x "))



# ==============================================================================
# 6. Feature filtering
# ==============================================================================
message("Filtering features")
taxa_before <- ncol(taxa_mat)
pathways_before_clean_filtering <- ncol(pathway_mat)

message("Removing special and taxon-stratified HUMAnN pathway features")
pathway_mat <- filter_clean_pathway_matrix(pathway_mat)
check_clean_pathway_matrix(pathway_mat)

#Remove samples with zero signal
taxa_mat <- remove_zero_sum_features(taxa_mat)
pathway_mat <- remove_zero_sum_features(pathway_mat)

# Keep only those taxes/pathways present in at least min_prevalence (10%) samples
taxa_mat <- filter_by_prevalence(
  mat = taxa_mat,
  min_prevalence = min_prevalence
)

pathway_mat <- filter_by_prevalence(
  mat = pathway_mat,
  min_prevalence = min_prevalence
)

taxa_after <- ncol(taxa_mat)
pathways_after <- ncol(pathway_mat)

# Save filtered but untransformed matrices for alpha diversity analysis.
taxa_raw_filtered <- taxa_mat
pathway_raw_filtered <- pathway_mat

message("Taxa features before filtering: ", taxa_before)
message("Taxa features after filtering: ", taxa_after)
message("Pathway features before filtering: ", pathways_before)
message("Pathway features after filtering: ", pathways_after)

if (taxa_after == 0) {
  stop("No taxonomic features left after filtering.")
}

if (pathways_after == 0) {
  stop("No pathway features left after filtering.")
}


# ==============================================================================
# 7. Transformation
# ==============================================================================

message("Transforming abundance matrices...")

taxa_mat <- transform_abundance(
  mat = taxa_mat,
  method = transformation,
  pseudocount = pseudocount
)

pathway_mat <- transform_abundance(
  mat = pathway_mat,
  method = transformation,
  pseudocount = pseudocount
)


# ==============================================================================
# 8. Clean feature names
# ==============================================================================

message("Cleaning feature names...")

safe_taxa_names <- make_safe_feature_names(colnames(taxa_mat))
safe_pathway_names <- make_safe_feature_names(colnames(pathway_mat))

colnames(taxa_mat) <- make_safe_feature_names(colnames(taxa_mat))
colnames(pathway_mat) <- make_safe_feature_names(colnames(pathway_mat))

colnames(taxa_raw_filtered) <- safe_taxa_names
colnames(pathway_raw_filtered) <- safe_pathway_names
# ==============================================================================
# 9. Build preprocessing summary
# ==============================================================================

message("Building preprocessing summary...")

preprocessing_summary <- build_preprocessing_summary(
  metadata = metadata,
  taxa_before = taxa_before,
  taxa_after = taxa_after,
  pathways_before = pathways_before,
  pathways_after = pathways_after,
  min_prevalence = min_prevalence,
  transformation = transformation
)

preprocessing_summary <- preprocessing_summary %>%
  bind_rows(
    tibble(
      metric = c(
        "n_pathways_before_clean_filtering",
        "n_pathways_after_clean_filtering",
        "n_pathways_removed_by_clean_filtering"
      ),
      value = c(
        as.character(pathways_before_clean_filtering),
        as.character(pathways_after_clean_filtering),
        as.character(pathways_before_clean_filtering - pathways_after_clean_filtering)
      )
    )
  )

# ==============================================================================
# 10. Save outputs
# ==============================================================================

message("Saving processed data...")

save_matrix_csv(taxa_mat, output_files$taxa)
save_matrix_csv(pathway_mat, output_files$pathways)

save_matrix_csv(taxa_raw_filtered, output_files$taxa_raw_filtered)
save_matrix_csv(pathway_raw_filtered, output_files$pathways_raw_filtered)

readr::write_csv(metadata, output_files$metadata)
readr::write_csv(preprocessing_summary, output_files$summary)


message("Saved transformed taxa matrix: ", output_files$taxa)
message("Saved transformed pathway matrix: ", output_files$pathways)
message("Saved raw filtered taxa matrix: ", output_files$taxa_raw_filtered)
message("Saved raw filtered pathway matrix: ", output_files$pathways_raw_filtered)
message("Saved processed metadata: ", output_files$metadata)
message("Saved preprocessing summary: ", output_files$summary)

message("Done.")
sessionInfo()