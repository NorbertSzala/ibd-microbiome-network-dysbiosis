suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(tibble)
})


# Function: load_abundance_matrix
# ------------------------------------------------------------------------------
# Description:
#   Loads a processed abundance matrix from CSV. The first column is named `sample_id` and  other columsn contain features
#
# Arguments:
#   path:
#     Path to the processed abundance df
#
#   Numeric matrix with sample IDs as row names and features as columns.
# ------------------------------------------------------------------------------
load_abundance_matrix <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  if (!"sample_id" %in% colnames(df)) {
    stop("Abundance matrix must contain a `sample_id` column: ", path)
  }

  mat <- df %>%
    column_to_rownames("sample_id") %>%
    as.matrix()

  storage.mode(mat) <- "numeric"

  mat
}

# ------------------------------------------------------------------------------
# Function: load_abundance_table
# ------------------------------------------------------------------------------
# Description:
#   Loads a processed abundance table from CSV. The first column must be named
#   `sample_id`, and all remaining columns are interpreted as abundance features.
#
# Arguments:
#   path:
#     Path to the processed abundance CSV file.
#
# Returns:
#   Tibble with `sample_id` column and abundance feature columns.
# ------------------------------------------------------------------------------
load_abundance_table <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  if (!"sample_id" %in% colnames(df)) {
    stop("Abundance table must contain a `sample_id` column: ", path)
  }

  df %>%
    dplyr::mutate(sample_id = as.character(sample_id))
}


# ------------------------------------------------------------------------------
# Function: load_metadata
# ------------------------------------------------------------------------------
# Description: 
#   Load preprocessed CSV metadata and checks whether required columns are present
#
# Arguments:
#   path:
#     Path to the processed csv data
#
# Returns:
#  A tibble containing metadata
# ------------------------------------------------------------------------------
load_metadata <- function(path) {
  metadata <- readr::read_csv(path, show_col_types=FALSE)

  required_columns <- c("sample_id", "disease_status")
  missing_columns <- setdiff(required_columns, colnames(metadata))

  if (length(missing_columns) > 0){
    stop("Metadata is missing required columns:", paste(missing_columns, collapse = ", "))
  } 

  metadata %>% 
    mutate(
        sample_id = as.character(sample_id),
        disease_status = factor(disease_status, levels = c("healthy", "IBD"))
    )
}


# ------------------------------------------------------------------------------
# Function: match_samples
# ------------------------------------------------------------------------------
# Description:
#   Keeps only samples present in both abundance matrix and metadata, and orders matrix rows to match metadata.
#
# Arguments:
#   mat:
#     Numeric sample by feature matrix. (innerjoin does not work on matrix)
#
#   metadata:
#     Metadata tibble containing `sample_id`.
#
# Returns:
#   List with matched matrix and metadata.
# ------------------------------------------------------------------------------
match_samples <- function(mat, metadata) {
    common_samples <- intersect(rownames(mat), metadata$sample_id)

    if (length(common_samples)==0) {
        stop("No common samples found between abundance matrix and metadata")
    }

    metadata_matched <- metadata %>% 
        filter(sample_id %in% common_samples) %>% 
        arrange(match(sample_id, common_samples))
    
    mat_matched <- mat[metadata_matched$sample_id, , drop=FALSE]

    list(mat=mat_matched, 
        metadata = metadata_matched)
}



# ------------------------------------------------------------------------------
# Function: clean_taxa_label
# ------------------------------------------------------------------------------
clean_taxa_label <- function(feature, max_words = 3, max_width = 45) {
  cleaned <- feature %>%
    stringr::str_replace("^_+", "") %>%
    stringr::str_replace_all("__+", " ") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish()

  shortened <- purrr::map_chr(
    stringr::str_split(cleaned, "\\s+"),
    ~ paste(head(.x, max_words), collapse = " ")
  )

  stringr::str_trunc(shortened, width = max_width)
}
# ------------------------------------------------------------------------------
# Function: clean pathway label
# ------------------------------------------------------------------------------
clean_pathway_label <- function(feature, max_width = 65) {
  feature %>%
    stringr::str_replace("(\\||_|__)(unclassified)$", "") %>%
    stringr::str_replace("^_+", "") %>%
    stringr::str_replace_all("__+", " ") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_trunc(width = max_width)
}

make_plot_label <- function(feature, feature_set) {
  if (feature_set == "taxa") {
    clean_taxa_label(feature)
  } else {
    clean_pathway_label(feature)
  }
}