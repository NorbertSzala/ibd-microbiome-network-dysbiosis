suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(tidyverse)
  library(readr)
  library(tibble)
  library(stringr)
})

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
extract_abundance_matrix <- function(se, assay_name) {
  assay_names <- SummarizedExperiment::assayNames(se)

  if (!assay_name %in% assay_names) {
    stop(
      "Assay `", assay_name, "` not found. Available assays: ",
      paste(assay_names, collapse = ", ")
    )
  }

  mat <- SummarizedExperiment::assay(se, assay_name)
  mat <- as.matrix(mat)
  mat <- t(mat) #TRanspose to get samples as rows and features as columns

  storage.mode(mat) <- "numeric" #change data type in matrix

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
  
  if (!sample_id_column %in% colnames(metadata)) {
    stop("Metadata must contain sample ID column: ", sample_id_column)
  }
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
# Function: match_taxa_pathway_samples
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
match_taxa_pathway_samples <- function(taxa_mat, pathway_mat, metadata) {
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

  # Make sure the matrix rows are in the same order as metadata rows
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
#     - "log1p": log(1 + x). --- log() is not possible because a lot of values in the matrices are zero.
#
# Returns:
#   Transformed numeric matrix.
# ------------------------------------------------------------------------------
transform_abundance <- function(mat, method = "log1p") {
  method <- as.character(method)

  if (method == "log1p") {
    return(log1p(mat))
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
# '__' Splits different levels of details in the features
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


# ------------------------------------------------------------------------------
# Function: filter_clean_pathway_matrix
# ------------------------------------------------------------------------------
# Description:
#   Keeps only unstratified community-level HUMAnN pathway features.
#   Removes:
#     - special HUMAnN rows such as UNMAPPED, UNINTEGRATED, UNGROUPED,
#     - taxon-stratified pathway features containing taxonomic identifiers.
#
#   This prevents mixing community-level functional profiles with taxon-specific
#   pathway contributions.
#
# Arguments:
#   mat:
#     Sample-by-feature pathway abundance matrix.
#
# Returns:
#   Matrix containing only candidate community-level pathway features.
# ------------------------------------------------------------------------------
filter_clean_pathway_matrix <- function(mat) {
  feature_names <- colnames(mat)

  is_special <- stringr::str_detect(
    feature_names,
    "^(UNMAPPED|UNINTEGRATED|UNGROUPED|UNCLASSIFIED|UNKNOWN|NA)(_|$)"
  )

  is_taxon_stratified <- stringr::str_detect(
    feature_names,
    "(\\||_)(k|p|c|o|f|g|s|t)__"
  )

  is_unclassified_stratified <- stringr::str_detect(
    feature_names,
    "(\\||_|__)(unclassified)$"
  )

  remove_features <- is_special | is_taxon_stratified | is_unclassified_stratified
  keep_features <- feature_names[!remove_features]

  message("Pathway features before clean filtering: ", length(feature_names))
  message("Removed special HUMAnN features: ", sum(is_special))
  message("Removed taxon-stratified features: ", sum(is_taxon_stratified & !is_special))
  message("Removed unclassified-stratified features: ", sum(is_unclassified_stratified & !is_special & !is_taxon_stratified))
  message("Removed total: ", sum(remove_features))
  message("Pathway features after clean filtering: ", length(keep_features))

  if (length(keep_features) == 0) {
    stop("No pathway features left after filtering.")
  }

  mat[, keep_features, drop = FALSE]
}