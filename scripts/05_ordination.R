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

set.seed(123)



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


# ==============================================================================
# 3. Helper functions
# ==============================================================================

# Function: load_abundance_matrix
# ------------------------------------------------------------------------------
# Description:
#   Loads a processed abundance matrix from CSV. The first column is named `sample_id` and  other columsn contain features
#
# Arguments:
#   path:
#     Path to the processed abundance matrix.
#
# Returns:
#   A data frame with `sample_id` and feature columns.
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
#     Numeric sample by feature matrix.
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
# Function: remove_zero_variances
# ------------------------------------------------------------------------------
# Description:
#   Removes features with zero variances across samples (those features that does not change between samples). Those features doeas contribute to PCA and may cause some issues
#
# Arguments:
#   mat:
#       numeric sample by feature matrix
#
# Returns:
#   Matrix without zer-variance columns
# ------------------------------------------------------------------------------
remove_zero_variances<- function(mat) {
    variances <- apply(mat, 2, var, na.rm=  TRUE) #Calculate variance for every column (feature). 
    keep <- variances > 0 & !is.na(variances)

    mat[, keep, drop = FALSE]
}


# ------------------------------------------------------------------------------
# Function: remove_empty_samples
# ------------------------------------------------------------------------------
# Description:
#   Removes samples with total abundance equal to zero. Those samples cant be used with bray curtis distance
# Arguments:
#   mat:
#     Numeric sample by feature matrix.
#
#   metadata:
#     Metadata table matched to matrix rows.
#
# Returns:
#   List with filtered matrix and metadata.
# ------------------------------------------------------------------------------
remove_empty_samples <- function(mat, metadata) {
  row_sums <- rowSums(mat, na.rm = TRUE)
  keep <- row_sums > 0 & !is.na(row_sums)

  n_removed <- sum(!keep)

  if (n_removed > 0) {
    message("Removing empty samples before Bray Curtis analysis: ", n_removed)
  }

  mat_filtered <- mat[keep, , drop = FALSE]
  metadata_filtered <- metadata[keep, , drop = FALSE]

  if (nrow(mat_filtered) == 0) {
    stop("All samples were removed as empty. Cannot continue ordination.")
  }

  list(
    mat = mat_filtered,
    metadata = metadata_filtered
  )
}

# ------------------------------------------------------------------------------
# Function: run_pca
# ------------------------------------------------------------------------------
# Description: 
#   Runs PCA on transformet abundance matrix and returns sample scores and variance explained by the first component
#
# Arguments:
#   mat:
#       numeric sample by feature matrix
#   metadata:
#       Metadata tibble containing 'sample id' and 'disease_status"
#   feature_set:
#       taxa/pathwyas for example
#
# Returns:
#  List with PCA scores and variance table
# ------------------------------------------------------------------------------
run_pca <- function(mat, metadata, feature_set) {
    mat <- remove_zero_variances(mat)
    
    filtered <- remove_empty_samples(mat, metadata)
    mat <- filtered$mat
    metadata <- filtered$metadata
    
    if (any(is.na(mat))) {
      stop("Missing values detected in matrix before PCA.")
    }

    pca <- prcomp(mat, center = TRUE, scale. = FALSE)

    explained <- (pca$sdev^2)/sum(pca$sdev^2) #component variance = pca$sdev^2

    scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>% 
      rownames_to_column("sample_id") %>% 
      left_join(metadata, by = "sample_id")%>% 
      mutate(feature_set = feature_set)

    variance <- tibble(
      feature_set = feature_set, 
      method = "PCA",
      axis = paste0("PC", seq_along(explained)),
      variance_explained = explained
    )

  list(scores = scores,
    variance = variance)
}


# ------------------------------------------------------------------------------
# Function: plot_pca
# ------------------------------------------------------------------------------
# Description:
#   Creates PCA scatter plot colored by disease status.
#
# Arguments:
#   pca_result:
#     Output list returned by run_pca() function.
#
#   title:
#     Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_pca <- function(pca_result, title) {
  var_pc1 <- pca_result$variance %>% 
    filter(axis == "PC1") %>% 
    pull(variance_explained) %>% 
    first()

    var_pc2 <- pca_result$variance %>% 
      filter(axis == "PC2") %>% 
      pull(variance_explained) %>% 
      first()
    
    ggplot(
      pca_result$scores,
      aes(x = PC1, y=PC2, color = disease_status)
    ) + 
    geom_point(alpha = 0.75, size = 2) +
    theme_minimal(base_size = 12) +
    labs(
      title = title,
      x = paste0("PC1 (", round(var_pc1*100, 1), "%)"),
      y = paste0("PC2 (", round(var_pc2*100, 1), "%)",
      color = "Disease status"
    ))
}


# ------------------------------------------------------------------------------
# Function: run_pcoa_bray
# ------------------------------------------------------------------------------
# Description:
#   Runs PCoA on Bray-Curtis distance matrix and returns sample coordinates and variance explained by ordination axes.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   metadata:
#     Metadata tibble containing `sample_id` and `disease_status`.
#
#   feature_set:
#       taxa/pathwyas for example
#
# Returns:
#  List with PCA scores and variance table

# We use Bray-Curtis as metric to count distances between each pair of samples. Then PCoA is made on this new matrix
# PCA = matrix (samples x feature)
# PCoA = matrix (samples x samples)
# ------------------------------------------------------------------------------
run_pcoa_bray <- function(mat, metadata, feature_set) {
  mat <- remove_zero_variances(mat)

  filtered <- remove_empty_samples(mat, metadata)
  mat <- filtered$mat
  metadata <- filtered$metadata


  if (any(mat < 0, na.rm = TRUE)) {
    stop(
      "Negative values detected. Bray-Curtis distance requires non-negative data. ",
      "Use log1p-transformed data, not CLR-transformed data, for Bray-Curtis PCoA."
    )
  }



  bray_dist <- vegan::vegdist(mat, method = "bray")

  if (any(is.na(bray_dist))) {
    stop(
      "Bray-Curtis distance contains NA values after removing empty samples. ",
      "Check whether the matrix contains missing values or invalid rows."
    )
  }

  pcoa <- cmdscale(
    bray_dist,
    eig = TRUE,
    k = 2
  )

  scores <- as.data.frame(pcoa$points) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("sample_id") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(feature_set = feature_set)

  positive_eig <- pcoa$eig[pcoa$eig > 0]
  variance_explained <- positive_eig / sum(positive_eig)

  variance <- tibble(
    feature_set = feature_set,
    method = "PCoA_BrayCurtis",
    axis = paste0("PCoA", seq_along(variance_explained)),
    variance_explained = variance_explained
  )

  list(
    scores = scores,
    variance = variance
  )
}


# ------------------------------------------------------------------------------
# Function: plot_pcoa
# ------------------------------------------------------------------------------
# Description:
#   Creates PCoA scatter plot colored by disease status.
#
# Arguments:
#   pcoa_result:
#     Output list returned by `run_pcoa_bray()`.
#
#   title:
#     Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_pcoa <- function(pcoa_result, title) {
  var_axis1 <- pcoa_result$variance %>%
    filter(axis == "PCoA1") %>%
    pull(variance_explained) %>%
    first()

  var_axis2 <- pcoa_result$variance %>%
    filter(axis == "PCoA2") %>%
    pull(variance_explained) %>%
    first()

  ggplot(
    pcoa_result$scores,
    aes(x = PCoA1, y = PCoA2, color = disease_status)
  ) +
    geom_point(alpha = 0.75, size = 2) +
    theme_minimal(base_size = 12) +
    labs(
      title = title,
      x = paste0("PCoA1 (", round(var_axis1 *100, 1), "%)"),
      y = paste0("PCoA2 (", round(var_axis2 *100, 1), "%)"),
      color = "Disease status"
    )
}


# ------------------------------------------------------------------------------
# Function: run_permanova
# ------------------------------------------------------------------------------
# Description:
#   Runs PERMANOVA using Bray-Curtis distance to test whether disease status explains global differences in microbiome composition.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   metadata:
#     Metadata tibble containing `disease_status`.
#
#   feature_set:
#     taxa/pathways
#
# Returns:
#   Tidy tibble with PERMANOVA results.
# ------------------------------------------------------------------------------
run_permanova <- function(mat, metadata, feature_set) {
  mat <- remove_zero_variances(mat)

  filtered <- remove_empty_samples(mat, metadata)
  mat <- filtered$mat
  metadata <- filtered$metadata

  if (any(mat < 0, na.rm = TRUE)) {
    stop(
      "Negative values detected. Bray-Curtis PERMANOVA requires non-negative data. ",
      "Use log1p-transformed data, not CLR-transformed data."
    )
  }

  bray_dist <- vegan::vegdist(mat, method = "bray")

  if (any(is.na(bray_dist))) {
    stop(
      "Bray-Curtis distance contains NA values after removing empty samples. ",
      "Check whether the matrix contains missing values or invalid rows."
    )
  }
  
  permanova <- vegan::adonis2(
    bray_dist ~ disease_status,
    data = metadata,
    permutations = 999
  )

  as.data.frame(permanova) %>%
    rownames_to_column("term") %>%
    mutate(
      feature_set = feature_set,
      distance = "Bray-Curtis"
    ) %>%
    select(feature_set, distance, term, Df, SumOfSqs, R2, F, `Pr(>F)`) %>%
    rename(
      p_value = `Pr(>F)`
    )
}




# ==============================================================================
# 4. Load data
# ==============================================================================
message("Loading input data")
taxa_mat <- load_abundance_matrix(input_files$taxa)
pathway_mat <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)


# ==============================================================================
# 5. Match samples
# ==============================================================================
message("Matching samples")
taxa_matched <- match_samples(taxa_mat, metadata)
pathway_matched <- match_samples(pathway_mat, metadata)

taxa_mat <- taxa_matched$mat
taxa_metadata <- taxa_matched$metadata

pathway_mat <- pathway_matched$mat
pathway_metadata <- pathway_matched$metadata


# ==============================================================================
# 6. PCA
# ==============================================================================
message("Running PCA")

pca_taxa <- run_pca(
  mat = taxa_mat,
  metadata = taxa_metadata,
  feature_set = "taxa"
)

pca_pathways <- run_pca(
  mat = pathway_mat,
  metadata = pathway_metadata,
  feature_set = "pathways"
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
# 7. PCoA
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
# 8. PERMANOVA
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
# 9. Save outputs
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
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = output_files$pca_pathways,
  plot = pca_pathway_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = output_files$pcoa_taxa,
  plot = pcoa_taxa_plot,
  width = 6,
  height = 5,
  dpi = 300
)

ggsave(
  filename = output_files$pcoa_pathways,
  plot = pcoa_pathway_plot,
  width = 6,
  height = 5,
  dpi = 300
)

message("Done")
sessionInfo()