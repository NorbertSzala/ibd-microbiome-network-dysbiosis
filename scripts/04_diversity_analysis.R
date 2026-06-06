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

set.seed(123)


# ==============================================================================
# 2. Paths
# ==============================================================================

# When running with Snakemake, paths can be passed automatically.
# For manual runs, define paths here.

if (exists("snakemake")) {
  input_files <- snakemake@input
  output_files <- snakemake@output
} else {
    stop("This script should be run through Snakemake")
}

walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. Helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
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
    readr::read_csv(path, show_col_types=FALSE)
}


# ------------------------------------------------------------------------------
# Function: load_metadata
# ------------------------------------------------------------------------------
# Description: 
#   Load preprocessed CSV metadata and checks whether required columns are present
#   #
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

  metadata
}



# ------------------------------------------------------------------------------
# Function:  abundance_df_to_matrix
# ------------------------------------------------------------------------------
# Description:
#   Converts processed abundance df into a numeric matrix with samples as rows and features as columns (transpose)
#
# Arguments:
#   abundance_df:
#     Data frame with 'sample_id' column and feature columns
#
# Returns:
#  Numeric matrix with rownames corresponding t sample IDs
# ------------------------------------------------------------------------------
abundance_df_to_matrix<- function(abundance_df) {
  if (!"sample_id" %in% colnames(abundance_df)) {
    stop("Abundance column must contain 'sample_id' in columns")
  }
  mat <- abundance_df %>% 
    column_to_rownames("sample_id") %>% 
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat
}


# ------------------------------------------------------------------------------
# Function: calculate_diversity_metrics_vegan
# ------------------------------------------------------------------------------
# Description:
#   Calculates alpha diversity metrics using the vegan package: Richness (specnumber), Shannon/Simpson/inverse Simpson (diversity), Chao1 (estimateR - on raw filtered data).
#
# Arguments:
#   mat:
#     Numeric sample-by-feature abundance matrix.
#
#   feature_set:
#     Name of the feature set, for example "taxa" or "pathways".
#
#   calculate_chao1:
#     Logical value. Chao1 should be calculated only on non-transformed abundance/count data.
#
# Returns:
#   Long-format tibble with columns: sample_id, feature_set, metric, value.
# ------------------------------------------------------------------------------
calculate_diversity_metrics_vegan<- function(mat, feature_set, calculate_chao1=FALSE) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  if (any(mat<0, na.rm = TRUE)) {
    stop(
      "Negative values detected in abundance matrix. Diversity metrics should be calculated only on row filtered non negative data"
    )
  }

  diversity_wide <- tibble(
    sample_id = rownames(mat),
    richness  = vegan::specnumber(mat), # typical richness, but counts even rarely 'organisms' or artefacts
    richness_threshold_1e5 = rowSums(mat > 1e-5, na.rm = TRUE),
    shannon = vegan::diversity(mat, index=  'shannon'),
    simpson = vegan::diversity(mat, index = "simpson"),
    inverse_simpson = vegan::diversity(mat, index = "invsimpson")
  ) %>% 
    mutate(
      hill_q0 = richness,
      hill_q1 = exp(shannon),
      hill_q2 = inverse_simpson, 
        pielou_evenness = ifelse(
          richness > 1,
          shannon / log(richness),
          NA_real_ 
        )
    )

  #Impossible when data is normalized
  if (calculate_chao1){
    chao_estimates <- vegan::estimateR(t(mat))
    chao1 <- chao_estimates["S.chao1", ]
    chao1_se <- chao_estimates['se.chao1', ]
    
    diversity_wide <- diversity_wide %>% 
      mutate(chao1 = as.numeric(chao1[.data$sample_id]),
            chao1_se = as.numeric(chao1_se[.data$sample_id]))
  }

  diversity_wide %>% 
    pivot_longer(
      cols = -sample_id,
      names_to = "metric",
      values_to = "value"
      ) %>% 
      mutate(feature_set = feature_set)
}



# ------------------------------------------------------------------------------
# Function: join_metadata
# ------------------------------------------------------------------------------
# Description:
#   Adds disease status metadata to diversity results and ensures that only samples with available metadata are retained
#   
#
# Arguments:
#   diversity_df:
#     long format table
#   metadata:
#     metadata table containing 'sample_id' and 'disease_status'
# Returns:
#  Diversity table with disease status annotation
# ------------------------------------------------------------------------------
join_metadata <- function(diversity_df, metadata) {
  diversity_df %>% 
    left_join(
      metadata %>% select(sample_id, disease_status),
      by="sample_id"
    ) %>% 
    filter(!is.na(disease_status))
}


# ------------------------------------------------------------------------------
# Function: run_wilcoxon_tests
# ------------------------------------------------------------------------------
# Description: 
#   Runs wilcoxon rank-sum tests comparing diversity values between healthy and IBD samples for each feature set and diversity metrics
#
# Arguments:
#   diversity_df:
#     long-format diversity table with columns: feature_set, metric, value, disease_status
#
# Returns:
#  tibble with test statistics, raw p-values and FDR adjusted p-values
# ------------------------------------------------------------------------------
run_wilcoxon_tests <- function(diversity_df) {
  diversity_df %>% 
    group_by(feature_set, metric) %>% 
    summarize(
      n_healthy = sum(disease_status =='healthy'),
      n_IBD = sum(disease_status == 'IBD'),
      median_healthy = median(value[disease_status == 'healthy'], na.rm = TRUE),
      median_IBD = median(value[disease_status == 'IBD'], na.rm = TRUE),
      difference_median = median_IBD - median_healthy,
      p_value = wilcox.test(
        value ~ disease_status,
        data = pick(everything()), 
        exact = FALSE)$p.value,
        .groups = "drop"
      ) %>% mutate (
        p_adj = p.adjust(p_value, method = "BH")
      ) %>% 
      arrange(feature_set, p_adj)
}


# ------------------------------------------------------------------------------
# Function: create_diversity_plot
# ------------------------------------------------------------------------------
# Description:
#   Creates boxplots with jittered points for diversity metrics, comparing healthy and IBD samples.
#
# Arguments:
#   diversity_df:
#     Long-format diversity table with disease status.
#
#   selected_feature_set:
#     Feature set to plot, for example "taxa" or "pathways".
#
#   title:
#     Plot title.
#
# Returns:
#   A ggplot object.
# ------------------------------------------------------------------------------
create_diversity_plot <- function(diversity_df, selected_feature_set, title) {
  diversity_df %>%
    filter(feature_set == selected_feature_set, is.finite(value)) %>%
    mutate(
      disease_status = factor(disease_status, levels = c("healthy", "IBD")),
      metric = factor(
        metric,
        levels = c(
          "richness", "richness_threshold_1e5", "hill_q0", "hill_q1", "hill_q2", "shannon", "simpson","inverse_simpson", "pielou_evenness"
        )
      )
    ) %>%
    ggplot(aes(x = disease_status, y = value, fill = disease_status)) +
    geom_violin(trim = FALSE, alpha = 0.25, linewidth = 0.3, color = "grey40") +
    geom_boxplot(outlier.shape = NA, linewidth = 0.4, width = 0.3, alpha = 0.7) +
    geom_jitter(aes(color = disease_status), width = 0.05, alpha = 0.25, size = 0.6, show.legend=  FALSE) +
    facet_wrap(~ metric, scales = "free_y", ncol=  3) +
    scale_fill_manual(values = c("healthy" = "#4C72B0", "IBD" = "#C44E52")) +
    scale_color_manual(values = c("healthy" = "#4C72B0", "IBD" = "#C44E52")) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = 'none', panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(), strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold")
    )+
    labs(
      title = title,
      x = NULL,
      y = "Alpha diversity value"
    )
}




# ==============================================================================
# 4. Load data
# ==============================================================================

message("Loading input data...")
taxa <- load_abundance_matrix(input_files$taxa)
pathways <- load_abundance_matrix(input_files$pathways)
metadata <- load_metadata(input_files$metadata)

# ==============================================================================
# 5. Prepare matrices
# ==============================================================================

message("Preparing analysis")
taxa_mat <- abundance_df_to_matrix(taxa)
pathways_mat <- abundance_df_to_matrix(pathways)


# ==============================================================================
# 6. Calculate diversity metrics
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
# 7. Statistical tests
# ==============================================================================
message("Running Wilcoxon tests")
diversity_tests <- run_wilcoxon_tests(diversity_results)


# ==============================================================================
# 8. Create plots
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
# 9. Save outputs
# ==============================================================================
message("Saving outputs...")

readr::write_csv(diversity_results, output_files$diversity_table)
readr::write_csv(diversity_tests, output_files$test_table)

ggsave(
  filename = output_files$taxa_plot,
  plot = taxa_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = output_files$pathway_plot,
  plot = pathways_plot,
  width = 8,
  height = 5,
  dpi = 300
)


message("Finished successfully.")
sessionInfo()