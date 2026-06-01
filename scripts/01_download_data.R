#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_download_data.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 01-06-2026
#
# Description:
#   First step of project analysis. Script donwloads the data from curatedMetagenomicData database and saves in proper location
#
# Inputs:
#   No input
#
# Outputs:
#    taxa = "results/raw/taxa.rds",
#    pathways = "results/raw/pathways.rds",
#    metadata = "results/raw/metadata.csv"
#
# Usage:
#   Rscript 01_download_data.R
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
  output_files <- snakemake@output
  params <- snakemake@params
} else {
  output_files <- list(
    taxa = "results/raw/taxa.rds",
    pathways = "results/raw/pathways.rds",
    metadata = "results/raw/metadata.csv"
  )

  params <- list(
    dataset = "LloydPrice_2019"
  )
}


# Create output directories if needed
walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. Download data
# ==============================================================================
message("Downloading data from curatedMetagenomicData...")






# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------

message("Saving raw objects...")