#!/usr/bin/env Rscript

# ==============================================================================
# Script: 00_template.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: <your_name>
# Date: <YYYY-MM-DD>
#
# Description:
#   Short description of what this script does.
#
# Inputs:
#   - input_file_1: description
#   - input_file_2: description
#
# Outputs:
#   - output_file_1: description
#   - output_file_2: description
#
# Usage:
#   Rscript scripts/00_template.R
#
# Notes:
#   This script is intended to be run directly or through Snakemake.
# ==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
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
  params <- snakemake@params
} else {
  input_files <- list(
    input_1 = "data/raw/input_file.csv"
  )

  output_files <- list(
    output_1 = "results/tables/output_file.csv"
  )

  params <- list(
    min_prevalence = 0.10
  )
}


# Create output directories if needed
walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. Helper functions
# ==============================================================================

check_file_exists <- function(path) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path)
  }
}


save_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, path)
  message("Saved: ", path)
}


filter_by_prevalence <- function(mat, min_prevalence = 0.10) {
  # mat: samples x features abundance matrix
  prevalence <- colMeans(mat > 0, na.rm = TRUE)
  mat[, prevalence >= min_prevalence, drop = FALSE]
}


log1p_transform <- function(mat) {
  log1p(mat)
}


clr_transform <- function(mat, pseudocount = 1e-6) {
  # Centered log-ratio transformation
  # Assumes samples are rows and features are columns

  mat <- as.matrix(mat)
  mat <- mat + pseudocount

  log_mat <- log(mat)
  gm <- rowMeans(log_mat)

  sweep(log_mat, 1, gm, FUN = "-")
}


# ==============================================================================
# 4. Load data
# ==============================================================================

message("Loading input data...")

# Example:
# check_file_exists(input_files$input_1)
# df <- read_csv(input_files$input_1, show_col_types = FALSE)


# ==============================================================================
# 5. Main analysis
# ==============================================================================

message("Running analysis...")

# Main code goes here.


# ==============================================================================
# 6. Save outputs
# ==============================================================================

message("Saving outputs...")

# Example:
# save_table(results_df, output_files$output_1)


# ==============================================================================
# 7. Session info
# ==============================================================================

message("Finished successfully.")
sessionInfo()