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

set.seed(123)


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
# 3. Helper functions
# ==============================================================================
# ------------------------------------------------------------------------------
# Function: check_file_exists
# ------------------------------------------------------------------------------
# Description:
#   Checks whether a required input file exists. If the file is missing, the function stops the script with an informative error message.

# Arguments:
#   path:
#     Path

# Returns:
#   Invisibly returns TRUE if the file exists.
# ------------------------------------------------------------------------------
check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file does not exist: ", path)
  }
  invisible(TRUE)
}


# ------------------------------------------------------------------------------
# Function: extract_abundance_matrix
# ------------------------------------------------------------------------------
# Description:
#   Extracts the first assay from a SummarizedExperiment object and converts it into a sample-by-feature matrix.
#
# Arguments:
#   se:
#     A SummarizedExperiment or TreeSummarizedExperiment object.
#
# Returns:
#   A numeric matrix with samples as rows and taxa/pathways as columns.
#------------------------------------------------------------------------------
extract_abundance_matrix <- function(se) {
  assay_names <- assayNames(se)

  if (length(assay_names) == 0) {
    stop("The SummarizedExperiment object contains no assays.")
  }

  mat <- assay(se, assay_names[1])

  mat <- as.matrix(mat)
  mat <- t(mat)

  storage.mode(mat) <- "numeric"

  mat
}


# ------------------------------------------------------------------------------
# Function: standardize_sample_id
# ------------------------------------------------------------------------------
# Description:
#   Cleans sample ids -  removes whitespaces
# ------------------------------------------------------------------------------
standardize_sample_id <- function(sample_ids) {
  sample_ids %>%
    as.character() %>%
    str_trim()
}


# ------------------------------------------------------------------------------
# Function: validate_metadata_columns
# ------------------------------------------------------------------------------
# Description:
#   Checks whether the metadata table contains user-specified (in config.yaml) columns required for preprocessing.
#
# Arguments:
#   metadata:
#     Data frame with sample metadata.
#
#   sample_id_column:
#     Name of the column containing sample identifiers.
#
#   disease_column:
#     Name of the column containing disease/status labels.
#
# Returns:
#   Invisibly returns TRUE if all required columns exist.
# ------------------------------------------------------------------------------
validate_metadata_columns <- function(metadata, sample_id_column, disease_column) {
  required_columns <- c(sample_id_column, disease_column)

  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0) {
    stop(
      "Missing required metadata columns: ",
      paste(missing_columns, collapse = ", "),
      "\nAvailable columns are: ",
      paste(colnames(metadata), collapse = ", ")
    )
  }

  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# Function: standardize_metadata_columns
# ------------------------------------------------------------------------------
# Description:
#   Renames user-specified metadata columns to standard internal names:  `sample_id` and `original_disease_label`
#
# Arguments:
#   metadata:
#     Data frame with sample metadata.
#
#   sample_id_column:
#     Name of the original metadata column containing sample identifiers.
#
#   disease_column:
#     Name of the original metadata column containing disease/status labels.
#
# Returns:
#   Metadata data frame with standardized `sample_id` and   original_disease_label` columns.
# ------------------------------------------------------------------------------
standardize_metadata_columns <- function(metadata, sample_id_column, disease_column) {
    metadata %>%
        mutate(
        sample_id = as.character(.data[[sample_id_column]]),
        original_disease_label = as.character(.data[[disease_column]])
    )
}


# ------------------------------------------------------------------------------
# Function: map_disease_status_explicit
# ------------------------------------------------------------------------------
# Description:
#   Maps disease labels to two analysis groups using labels provided in config.yaml. Samples matching `healthy_label` are assigned to "healthy", and samples matching `ibd_label` are assigned to "IBD". All other labels are set to NA and removed later.
#
# Arguments:
#   metadata:
#     Metadata data frame with `original_disease_label`.
#
#   healthy_label:
#     Label in metadata corresponding to healthy/control samples.
#
#   ibd_label:
#     Label in metadata corresponding to IBD samples.
#
# Returns:
#   Metadata data frame with a new factor column `disease_status`.
# ------------------------------------------------------------------------------
map_disease_status_explicit <- function(metadata, healthy_label, ibd_label) {
  healthy_label_clean <- str_to_lower(str_trim(healthy_label))
  ibd_label_clean <- str_to_lower(str_trim(ibd_label))

  metadata %>%
    mutate(
      disease_label_clean = str_to_lower(
        str_trim(.data$original_disease_label)
      ),
      disease_status = dplyr::case_when(
        .data$disease_label_clean == healthy_label_clean ~ "healthy",
        .data$disease_label_clean == ibd_label_clean ~ "IBD",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(.data$disease_status)) %>%
    mutate(
      disease_status = factor(.data$disease_status, levels = c("healthy", "IBD"))
    )
}


# ------------------------------------------------------------------------------
# Function: filter_metadata_to_groups
# ------------------------------------------------------------------------------
# Description:
#   Applies explicit metadata column selection and explicit disease label mapping.
# Arguments:
#   metadata:
#     Data frame with sample metadata.
#
#   sample_id_column:
#     Column name containing sample identifiers.
#
#   disease_column:
#     Column name containing disease/status labels.
#
#   healthy_label:
#     Label representing healthy/control samples.
#
#   ibd_label:
#     Label representing IBD samples.
#
# Returns:
#   Filtered metadata table containing only healthy and IBD samples.
# ------------------------------------------------------------------------------
filter_metadata_to_groups <- function(metadata, sample_id_column, disease_column, healthy_label, ibd_label
) {
  validate_metadata_columns(
    metadata = metadata,
    sample_id_column = sample_id_column,
    disease_column = disease_column
  )

  message("Using sample ID column: ", sample_id_column)
  message("Using disease/status column: ", disease_column)
  message("Healthy label: ", healthy_label)
  message("IBD label: ", ibd_label)

  metadata %>%
    standardize_metadata_columns(
      sample_id_column = sample_id_column,
      disease_column = disease_column
    ) %>%
    map_disease_status_explicit(
      healthy_label = healthy_label,
      ibd_label = ibd_label
    )
}


# ------------------------------------------------------------------------------
# Function: match_samples
# ------------------------------------------------------------------------------
# Description:
#   Keeps only samples present in all three objects: taxonomic matrix, pathway
#   matrix, and metadata. It also ensures that rows of both matrices and metadata are ordered in the same way.
#
# Arguments:
#   taxa_mat:
#     Sample-by-feature taxonomic matrix.
#
#   pathway_mat:
#     Sample-by-feature pathway matrix.
#
#   metadata:
#     Metadata data frame with a `sample_id` column.
#
# Returns:
#   A list containing matched `taxa`, `pathways`, and `metadata` objects.
# ------------------------------------------------------------------------------
match_samples <- function(taxa_mat, pathway_mat, metadata) {
  taxa_ids <- standardize_sample_id(rownames(taxa_mat))
  pathway_ids <- standardize_sample_id(rownames(pathway_mat))
  metadata_ids <- standardize_sample_id(metadata$sample_id)

  rownames(taxa_mat) <- taxa_ids
  rownames(pathway_mat) <- pathway_ids
  metadata$sample_id <- metadata_ids

  common_samples <- Reduce(
    intersect,
    list(
      rownames(taxa_mat),
      rownames(pathway_mat),
      metadata$sample_id
    )
  )

  if (length(common_samples) == 0) {
    stop("No common samples found between taxa, pathways, and metadata.")
  }

  message("Number of common samples: ", length(common_samples))

  taxa_mat <- taxa_mat[common_samples, , drop = FALSE]
  pathway_mat <- pathway_mat[common_samples, , drop = FALSE]

  metadata <- metadata %>%
    filter(.data$sample_id %in% common_samples) %>%
    arrange(match(.data$sample_id, common_samples))

  taxa_mat <- taxa_mat[metadata$sample_id, , drop = FALSE]
  pathway_mat <- pathway_mat[metadata$sample_id, , drop = FALSE]

  list(
    taxa = taxa_mat,
    pathways = pathway_mat,
    metadata = metadata
  )
}


# ------------------------------------------------------------------------------
# Function: remove_zero_sum_features
# ------------------------------------------------------------------------------
# Description:
#   Removes features whose total abundance across all samples is zero.
#
# Arguments:
#   mat:
#     Sample-by-feature abundance matrix.
#
# Returns:
#   Matrix without zero-sum columns.
# ------------------------------------------------------------------------------
remove_zero_sum_features <- function(mat) {
  keep <- colSums(mat, na.rm = TRUE) > 0
  mat[, keep, drop = FALSE]
}


# ------------------------------------------------------------------------------
# Function: filter_by_prevalence
# ------------------------------------------------------------------------------
# Description:
#   Removes rare features based on prevalence across samples.
#
# Arguments:
#   mat:
#     Sample-by-feature abundance matrix.
#
#   min_prevalence:
#     Minimum fraction of samples in which a feature must be present. Example: 0.10 means that a feature must be non-zero in at least 10% of samples.
#
# Returns:
#   Filtered sample-by-feature matrix.
# ------------------------------------------------------------------------------
filter_by_prevalence <- function(mat, min_prevalence = 0.10) {
  if (min_prevalence < 0 || min_prevalence > 1) {
    stop("min_prevalence must be between 0 and 1.")
  }

  prevalence <- colMeans(mat > 0, na.rm = TRUE)
  keep <- prevalence >= min_prevalence

  mat[, keep, drop = FALSE]
}


# ------------------------------------------------------------------------------
# Function: clr_transform
# ------------------------------------------------------------------------------
# Description:
#   Applies centered log-ratio transformation to compositional abundance data.
#
# Arguments:
#   mat:
#     Sample-by-feature abundance matrix.
#
#   pseudocount:
#     Small value added to all abundances to avoid log(0).
#
# Returns:
#   CLR-transformed numeric matrix.
#
# Notes:
#   CLR transformation is often used for compositional microbiome data. However,  the result depends on the pseudocount and feature filtering strategy.
# ------------------------------------------------------------------------------
clr_transform <- function(mat, pseudocount = 1e-6) {
  mat <- as.matrix(mat)
  mat <- mat + pseudocount

  log_mat <- log(mat)
  geometric_mean_log <- rowMeans(log_mat)

  sweep(log_mat, 1, geometric_mean_log, FUN = "-")
}


# ------------------------------------------------------------------------------
# Function: transform_abundance
# ------------------------------------------------------------------------------
# Description:
#   Applies the selected transformation method to an abundance matrix.
#
# Arguments:
#   mat:
#     Sample-by-feature abundance matrix.
#
#   method:
#     Transformation method. Supported values:
#     - "log1p": log(1 + x)
#     - "clr": centered log-ratio transformation
#
# Returns:
#   Transformed numeric matrix.
# ------------------------------------------------------------------------------
transform_abundance <- function(mat, method = "log1p") {
  method <- as.character(method)

  if (method == "log1p") {
    return(log1p(mat))
  }

  if (method == "clr") {
    return(clr_transform(mat))
  }

  stop("Unknown transformation method: ", method)
}


# ------------------------------------------------------------------------------
# Function: make_safe_feature_names
# ------------------------------------------------------------------------------
# Description:
#   Converts feature names into safe, unique column names suitable for CSV downstream analyss
#
# Arguments:
#   feature_names:
#     Character vector of original feature names.
#
# Returns:
#   Character vector of cleaned and unique feature names.
#
# ------------------------------------------------------------------------------
make_safe_feature_names <- function(feature_names) {
  feature_names %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9_\\.\\-]", "_") %>%
    make.unique(sep = "__")
}


# ------------------------------------------------------------------------------
# Function: save_matrix_csv
# ------------------------------------------------------------------------------
# Description:
#   Saves a sample-by-feature matrix to CSV with sample identifiers stored in the first column named `sample_id`.
#
# Arguments:
#   mat:
#     Sample-by-feature matrix.
#
#   path:
#     Output CSV path.
#
# Returns:
#   Invisibly returns the written path.
# ------------------------------------------------------------------------------
save_matrix_csv <- function(mat, path) {
  output_df <- as.data.frame(mat) %>%
    rownames_to_column("sample_id")

  readr::write_csv(output_df, path)

  invisible(path)
}


# ------------------------------------------------------------------------------
# Function: build_preprocessing_summary
# ------------------------------------------------------------------------------
# Description:
#   Builds a compact summary table describing how many samples and features were
#   retained after preprocessing.
#
# Arguments:
#   metadata:
#     Final metadata table.
#
#   taxa_before:
#     Number of taxonomic features before filtering.
#
#   taxa_after:
#     Number of taxonomic features after filtering.
#
#   pathways_before:
#     Number of pathway features before filtering.
#
#   pathways_after:
#     Number of pathway features after filtering.
#
#   min_prevalence:
#     Prevalence threshold used for filtering.
#
#   transformation:
#     Transformation method used.
#
# Returns:
#   A tibble with preprocessing summary statistics.
# ------------------------------------------------------------------------------
build_preprocessing_summary <- function(
  metadata,
  taxa_before,
  taxa_after,
  pathways_before,
  pathways_after,
  min_prevalence,
  transformation
) {
  group_counts <- metadata %>%
    count(.data$disease_status, name = "n_samples") %>%
    mutate(disease_status = as.character(.data$disease_status))

  tibble(
    metric = c(
      "n_samples_total",
      "n_taxa_before_filtering",
      "n_taxa_after_filtering",
      "n_pathways_before_filtering",
      "n_pathways_after_filtering",
      "min_prevalence",
      "transformation"
    ),
    value = c(
      as.character(nrow(metadata)),
      as.character(taxa_before),
      as.character(taxa_after),
      as.character(pathways_before),
      as.character(pathways_after),
      as.character(min_prevalence),
      transformation
    )
  ) %>%
    bind_rows(
      group_counts %>%
        transmute(
          metric = paste0("n_samples_", .data$disease_status),
          value = as.character(.data$n_samples)
        )
    )
}


# ==============================================================================
# 4. Load data
# ==============================================================================

message("Checking input files...")

check_file_exists(input_files$taxa)
check_file_exists(input_files$pathways)
check_file_exists(input_files$metadata)

message("Loading raw objects...")

taxa_se <- readRDS(input_files$taxa)
pathway_se <- readRDS(input_files$pathways)
metadata <- readr::read_csv(input_files$metadata, show_col_types = FALSE)

if (!"sample_id" %in% colnames(metadata)) {
  stop("Metadata must contain a `sample_id` column.")
}


# ==============================================================================
# 5. Extract abundance matrices
# ==============================================================================

message("Extracting abundance matrices...")

taxa_mat <- extract_abundance_matrix(taxa_se)
pathway_mat <- extract_abundance_matrix(pathway_se)

message("Raw taxa matrix dimensions: ", paste(dim(taxa_mat), collapse = " x "))
message("Raw pathway matrix dimensions: ", paste(dim(pathway_mat), collapse = " x "))


# ==============================================================================
# 6. Metadata filtering and sample matching
# ==============================================================================

message("Filtering metadata to healthy and IBD samples...")

metadata <- filter_metadata_to_groups(
  metadata = metadata,
  sample_id_column = sample_id_column,
  disease_column = disease_column,
  healthy_label = healthy_label,
  ibd_label = ibd_label
)

matched <- match_samples(
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
# 7. Feature filtering
# ==============================================================================

message("Filtering features...")

taxa_before <- ncol(taxa_mat)
pathways_before <- ncol(pathway_mat)

taxa_mat <- remove_zero_sum_features(taxa_mat)
pathway_mat <- remove_zero_sum_features(pathway_mat)

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
# 8. Transformation
# ==============================================================================

message("Transforming abundance matrices...")

taxa_mat <- transform_abundance(
  mat = taxa_mat,
  method = transformation
)

pathway_mat <- transform_abundance(
  mat = pathway_mat,
  method = transformation
)


# ==============================================================================
# 9. Clean feature names
# ==============================================================================

message("Cleaning feature names...")

safe_taxa_names <- make_safe_feature_names(colnames(taxa_mat))
safe_pathway_names <- make_safe_feature_names(colnames(pathway_mat))

colnames(taxa_mat) <- make_safe_feature_names(colnames(taxa_mat))
colnames(pathway_mat) <- make_safe_feature_names(colnames(pathway_mat))

colnames(taxa_raw_filtered) <- safe_taxa_names
colnames(pathway_raw_filtered) <- safe_pathway_names
# ==============================================================================
# 10. Build preprocessing summary
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


# ==============================================================================
# 11. Save outputs
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