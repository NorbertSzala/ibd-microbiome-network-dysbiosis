#!/usr/bin/env Rscript

# ==============================================================================
# Script: 01_download_data.R
# Project: Functional vs Taxonomic Dysbiosis in IBD Gut Microbiomes
# Author: Norbert Szala
# Date: 01-06-2026
#
# Description:
#   Downloads taxonomic relative abundance profiles and metabolic pathway
#   abundance profiles from curatedMetagenomicData. The script saves raw
#   SummarizedExperiment objects and sample metadata for downstream preprocessing. (source: https://bioconductor.org/packages/3.16/data/experiment/manuals/curatedMetagenomicData/man/curatedMetagenomicData.pdf)
#
# Outputs:
#   - taxa.rds
#   - pathways.rds
#   - metadata.csv
#   - sample_summary.csv
#==============================================================================

# NOTE:
# Raw data does not mean raw reads/FASTQ. It means raw for this pipeline. The Data is already preprocessed and it contains relative abundance od features (taxa/metanolic pathways) in the samples.
# 
# The bacterial, fungal, and archaeal taxonomic abundances for each sample were calculated with MetaPhlAn3, and metabolic functional potential was calculated with HUMAnN3
# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(tidyverse)
})

source("scripts/functions/01_download_data_functions.R")


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
disease_col <- params$disease_column #TODO: Jaka disease column

message("Selected study: ", study_name)
message("Selected body site: ", target_body_site)

# ==============================================================================
# 3. Download data
# ==============================================================================
message("Downloading taxonomic data from curatedMetagenomicData...")

taxa_se <- get_first_resource(
  pattern = paste0(study_name, ".*relative_abundance"),
  rownames = "short"
)
message("Downloading metabolic pathway data from curatedMetagenomicData...")
pathway_se <- get_first_resource(
  pattern = paste0(study_name, ".*pathway_abundance"),
  rownames = 'short'
)

taxa_samples <- colnames(taxa_se)
pathway_samples <- colnames(pathway_se)

message("Taxa samples: ", length(taxa_samples))
message("Pathway samples: ", length(pathway_samples))
message("Shared samples: ", length(intersect(taxa_samples, pathway_samples)))


# ------------------------------------------------------------------------------
# 4. Extract metadata
# ------------------------------------------------------------------------------
message("Extracting metadata")
metadata <- clean_metadata(se = taxa_se, target_body_site = target_body_site)

sample_summary <- summarize_metadata(metadata, disease_col)

message("Number of metadata samples after body site filtering: ", nrow(metadata))


message("Filtered metadata samples in taxa: ",
        sum(metadata$sample_id %in% colnames(taxa_se)))

message("Filtered metadata samples in pathways: ",
        sum(metadata$sample_id %in% colnames(pathway_se)))

print(sample_summary)
# ------------------------------------------------------------------------------
# 5. Save outputs
# ------------------------------------------------------------------------------
message("Saving raw .rds objects and metadata...")

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