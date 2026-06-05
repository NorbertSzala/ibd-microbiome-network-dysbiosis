#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_download_data.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 01-06-2026
#
# Description:
#   Downloads taxonomic relative abundance profiles and functional pathway
#   abundance profiles from curatedMetagenomicData. The script saves raw
#   SummarizedExperiment objects and sample metadata for downstream preprocessing. (source: https://bioconductor.org/packages/3.16/data/experiment/manuals/curatedMetagenomicData/man/curatedMetagenomicData.pdf)
#
# Inputs:
#   No file inputs.
#
# Outputs:
#   - taxa.rds
#   - pathways.rds
#   - metadata.csv
#   - sample_summary.csv
#
# Usage:
#   Rscript scripts/01_download_data.R
#
# Notes:
#   This script can be run directly or through Snakemake.
# ==============================================================================


# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(tidyverse)
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
  stop("R cannot find snakemake output and params")
}

walk(dirname(unlist(output_files)), dir.create, recursive = TRUE, showWarnings = FALSE)

study_name <- params$study_name
target_body_site <- params$body_site

message("Selected study: ", study_name)
message("Selected body site: ", target_body_site)


# ==============================================================================
# 3. Helper functions
# ==============================================================================
get_first_resource <- function(pattern, rownames = NULL) {
  message("Searching resource: ", pattern)

  available <- curatedMetagenomicData(pattern = pattern, dryrun = TRUE)

  if (length(available)==0) {stop("No resources found for pattern: ", pattern)}

  message("Available resources: ")
  print(available)
  message("Downloading the most recent matching resource")

  if (is.null(rownames)) {# Then default is "long"
    resource_list <- curatedMetagenomicData(pattern=pattern, dryrun=FALSE)
  } else {
    resource_list <- curatedMetagenomicData(pattern = pattern, dryrun=FALSE, rownames = rownames)
  }

  if (length(resource_list)) == 0 {
    stop("Downloading resources from curatedMetagenomicData returned no object for pattern: ", pattern)
  }

  # curatedMetagenomicData returns a list
  resource_list[[1]]
}


clean_metadata <- function(se, target_body_site = "stool") {
  metadata <- as.data.frame(colData(se)) %>%
    tibble::rownames_to_column("sample_id")

  # For now, keep all columns but make sure sampleid is explicit
  if ("body_site" %in% colnames(metadata)) {
    metadata <- metadata %>%
      filter(body_site == target_body_site)
  } else {
    warning("Column 'body_site' not found in metadata. No body site filtering applied.")
  }

  metadata
}

summarize_metadata <- function(metadata) {
  target_columns <- c("disease", "study_condition", "disease_status", "diagnosis")

  disease_col <- intersect(target_columns, colnames(metadata))

  if (length(disease_col) == 0) {
    return(tibble(variable = "n_samples", value=nrow(metadata)))
  }

  disease_col <- disease_col[[1]]

  metadata %>%
    count(.data[[disease_col]], name = "n") %>%
    rename(group=1) %>%
    mutate(variable=disease_col) %>%
    select(variable, group, n)
}



# ==============================================================================
# 4. Download data
# ==============================================================================
message("Downloading taxonomic data from curatedMetagenomicData...")
taxa_se <- get_first_resource(
    pattern = paste0(study_name, ".relative_abundance"), 
    rownames = "short")

message("Downloading taxonomic data from curatedMetagenomicData...")
pathway_se <- get_first_resource(
  pattern = paste0(study_name, ".pathway_abundance")
)



# ------------------------------------------------------------------------------
# 5. Extract metadata
# ------------------------------------------------------------------------------
message("Extracting metadata")
metadata <- clean_metadata(se = taxa_se, target_body_site = target_body_site)

sample_summary <- summarize_metadata(metadata)

message("Number of metadata samples after body site filtering: ", nrow(metadata))


# ------------------------------------------------------------------------------
# 6. Save outputs
# ------------------------------------------------------------------------------
message("Saving raw objects and metadata...")

saveRDS(taxa_se, output_files$taxa)
saveRDS(pathway_se, output_files$pathways)

readr::write_csv(metadata, output_files$metadata)
readr::write_csv(sample_summary, output_files$summary)

message("Saved taxa object: ", output_files$taxa)
message("Saved pathway object: ", output_files$pathways)
message("Saved metadata: ", output_files$metadata)
message("Saved sample summary: ", output_files$summary)

message("Done.")
sessionInfo()

# library(curatedMetagenomicData)

# curatedMetagenomicData("LloydPrice_2019", dryrun = TRUE)


# library(curatedMetagenomicData)
# library(dplyr)

# sampleMetadata %>%
#   filter(grepl("Lloyd|IBD|inflammatory|bowel", study_name, ignore.case = TRUE)) %>%
#   distinct(study_name)