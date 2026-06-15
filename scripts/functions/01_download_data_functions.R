suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(tidyverse)
})

set.seed(123)



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

  if (length(resource_list) == 0) {
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