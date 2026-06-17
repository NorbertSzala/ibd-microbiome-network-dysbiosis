suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(SummarizedExperiment)
  library(tidyverse)
})


# Finds resources matching a curatedMetagenomicData pttern and downloads the first matching object.
# The rownames argument controls feature naming

get_first_resource <- function(pattern, rownames = NULL) {
  message("Searching resource: ", pattern)

  # Checks in dryrun what data are available
  available <- curatedMetagenomicData(pattern = pattern, dryrun = TRUE)

  if (length(available)==0) {stop("No resources found for pattern: ", pattern)}

  if (length(available) > 1) {
    warning("Multiple resources found. Using first: ", available[[1]])
  }

  message("Available resources: ")
  print(available)
  message("Downloading the first matching resource")

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


# Choose only samples coming from one tissue
#SE = SummarizedExperiment (object)
clean_metadata <- function(se, target_body_site = "stool") {
  metadata <- as.data.frame(colData(se)) %>%
    tibble::rownames_to_column("sample_id")

  # For now, keep all columns but make sure sampleid is explicit. Usually stool
  if ("body_site" %in% colnames(metadata)) {
    metadata <- metadata %>%
      filter(body_site == target_body_site)
  } else {
    warning("Column 'body_site' not found in metadata. No body site filtering applied.")
  }

  metadata
}

# Create a simple count table for the disease/group variab
summarize_metadata <- function(metadata, disease_col) {
  if (!disease_col %in% colnames(metadata)) {
    stop("Disease column not found in metadata: ", disease_col)
  }

  summary <- metadata %>%
    count(.data[[disease_col]], name = "n")

  colnames(summary)[1] <- "group"

  summary %>%
    mutate(variable = disease_col) %>%
    select(variable, group, n)
}