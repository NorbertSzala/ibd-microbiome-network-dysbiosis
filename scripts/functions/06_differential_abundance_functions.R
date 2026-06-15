# ==============================================================================
# 1. Setup
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
  library(forcats)
})


# ==============================================================================
# 3. Helper functions
# ==============================================================================



# ------------------------------------------------------------------------------
# Function: matrix_to_long_df
# ------------------------------------------------------------------------------
# Description:
#   Converts a sample-by-feature abundance matrix to long format. The output
#   contains one row per sample-feature pair.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   feature_set:
#     Name of the feature set, for example "taxa" or "pathways".
#
# Returns:
#   Long-format tibble with columns: sample_id, feature, value, feature_set.
# ------------------------------------------------------------------------------
matrix_to_long_df <- function(mat, feature_set) {
  mat %>% 
    as.data.frame() %>% 
    rownames_to_column("sample_id") %>% 
    pivot_longer(
      cols = -sample_id,
      names_to = "feature",
      values_to = "value"
    ) %>% 
    mutate(
      sample_id = as.character(sample_id),
      feature_set = feature_set
    )
}


# ------------------------------------------------------------------------------
# Function: wilcoxon_test
# ------------------------------------------------------------------------------
# Description:
#   Runs Wilcoxon rank-sum test safely. If the feature has insufficient variation
#   or only one disease group is present, returns NA instead of stopping the script.
#
# Arguments:
#   value:
#     Numeric vector with feature abundance values.
#
#   group:
#     Factor or character vector with disease groups.
#
# Returns:
#   Numeric p-value or NA.
# ------------------------------------------------------------------------------
wilcoxon_test <- function(value, group) {
  df <- tibble(value=value, group = group) %>% 
    filter(!is.na(value), !is.na(group))

  if (n_distinct(df$group) < 2) {
    return(NA_real_)
  }

  if (n_distinct(df$value) < 2 ) {
    return(NA_real_)
  }

  wilcox.test(
    value~group,
    data = df,
    exact = FALSE
  )$p.value
}



# ------------------------------------------------------------------------------
# Function: run_differential_abundance
# ------------------------------------------------------------------------------
# Description:
#   Performs differential abundance analysis for each feature using Wilcoxon
#   rank-sum test. P-values are adjusted using Benjamini-Hochberg FDR correction.
#
# Arguments:
#   long_df:
#     Long-format abundance table with columns:
#     sample_id, feature, value, feature_set, disease_status.
#
# Returns:
#   Tibble with one row per feature and differential abundance statistics.
# ------------------------------------------------------------------------------
run_differential_abundance <- function(long_df, p_adjust_method = "BH") {
  required_columns <- c("sample_id", "feature", "value", "feature_set", "disease_status")
  missing_columns <- setdiff(required_columns, colnames(long_df))

  if (length(missing_columns) > 0 ) {
    stop("Differential abundance input is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  long_df %>% 
    group_by(feature_set, feature) %>% 
    summarize(
      n_healthy = sum(disease_status == "healthy"),
      n_ibd = sum(disease_status == "IBD"),

      mean_healthy = mean(value[disease_status=='healthy'], na.rm = TRUE),
      mean_ibd = mean(value[disease_status == 'IBD'], na.rm = TRUE),

      median_healthy = median(value[disease_status == "healthy"], na.rm = TRUE),
      median_ibd = median(value[disease_status == "IBD"], na.rm = TRUE),

      median_nonzero_healthy = ifelse(
        sum(disease_status == "healthy" & value > 0, na.rm = TRUE) > 0,
        median(value[disease_status == "healthy" & value > 0], na.rm = TRUE),
        NA_real_
      ),

      median_nonzero_ibd = ifelse(
        sum(disease_status == "IBD" & value > 0, na.rm = TRUE) > 0,
        median(value[disease_status == "IBD" & value > 0], na.rm = TRUE),
        NA_real_
      ),

      difference_median_nonzero = median_nonzero_ibd - median_nonzero_healthy,

      difference_mean = mean_ibd - mean_healthy,
      difference_median = median_ibd - median_healthy,

      log2_fold_change = log2(
        (mean_ibd + 1e-6) / (mean_healthy + 1e-6)
      ),
      prevalence_healthy = mean(value[disease_status == "healthy"] > 0, na.rm = TRUE),
      prevalence_ibd = mean(value[disease_status == "IBD"] > 0, na.rm = TRUE),
      difference_prevalence = prevalence_ibd - prevalence_healthy,

      p_value = wilcoxon_test(
        value = value,
        group = disease_status
      ),


      .groups = "drop"
    ) %>%
    group_by(feature_set) %>%
    mutate(
      p_adj = p.adjust(p_value, method = p_adjust_method),
      direction = case_when(
        difference_median > 0 ~ "higher_in_IBD",
        difference_median < 0 ~ "lower_in_IBD",
        TRUE ~ "no_median_difference"
      ),
      prevalence_direction = case_when(
      difference_prevalence > 0 ~ "more_prevalent_in_IBD",
      difference_prevalence < 0 ~ "less_prevalent_in_IBD",
      TRUE ~ "no_prevalence_difference"
    )) %>%
    ungroup() %>%
    arrange(feature_set, p_adj, desc(abs(difference_median)))
}



# ------------------------------------------------------------------------------
# Function: get_top_features
# ------------------------------------------------------------------------------
# Description:
#   Selects top differential features based on adjusted p-value and absolute
#   median difference.
#
# Arguments:
#   results_df:
#     Differential abundance results table.
#
#   ntop:
#     Number of top features to select per feature set.
#
# Returns:
#   Tibble with top features.
# ------------------------------------------------------------------------------
get_top_features <- function(results_df, ntop = 20) {
  results_df %>% 
    filter(!is.na(p_adj)) %>% 
    group_by(feature_set) %>% 
    arrange(
      p_adj,
      desc(abs(difference_prevalence)),
      desc(abs(difference_median))
    , .by_group = TRUE) %>% 
    slice_head(n=ntop) %>% 
    ungroup()
}



# ------------------------------------------------------------------------------
# Function: plot_top_features
# ------------------------------------------------------------------------------
# Description:
#   Creates a bar plot of top differentially abundant features, showing median
#   abundance difference between IBD and healthy samples.
#
# Arguments:
#   results_df:
#     Differential abundance results table.
#
#   feature_set:
#     Feature set to plot: "taxa" or "pathways".
#
#   ntop:
#     Number of top features to display.
#
#   title:
#     Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_top_features <- function(results_df, selected_feature_set, title, ntop = 20) {
  plotdf <- results_df %>%
    filter(feature_set == selected_feature_set )%>% 
    filter(!is.na(p_adj)) %>% 
    arrange(
      p_adj,
      desc(abs(difference_prevalence)),
      desc(abs(difference_median))
    ) %>% 
    slice_head(n = ntop) %>% 
    mutate(
      feature_short = stringr::str_replace_all(feature, "_", " "),
      feature_short = stringr::str_trunc(feature_short, width = 70),
      feature_short = forcats::fct_reorder(feature_short, difference_median)
    )

  ggplot(
    plotdf, 
    aes(x = feature_short, y = difference_median, fill = direction)
  ) +
    geom_col() + 
    coord_flip() +
    theme_minimal(base_size = 11)+
    labs(
      title = title,
      x = NULL, y = "Median difference: IBD - healthy", fill = "Direction"
    )
}


shorten_feature_name <- function(feature, n_words = 3) {
  feature %>%
    stringr::str_replace("__.*$", "") %>%
    stringr::str_replace_all("[_\\-]+", " ") %>%
    stringr::str_squish() %>%
    stringr::word(1, n_words)
}

plot_top_prevalence_features <- function(results_df, selected_feature_set, title, ntop = 20, max_label_words_taxa = 2, max_label_words_pathways = 3) {
  plot_df <- results_df %>%
    filter(feature_set == selected_feature_set) %>%
    filter(!is.na(p_adj)) %>%
    arrange(p_adj, desc(abs(difference_prevalence))) %>%
    slice_head(n = ntop) %>%
    mutate(
      feature_short = ifelse(
        selected_feature_set == "taxa",
        shorten_feature_name(feature, n_words = max_label_words_taxa),
        shorten_feature_name(feature, n_words = max_label_words_pathways)
      ),
      feature_short = forcats::fct_reorder(feature_short, difference_prevalence),
      prevalence_direction = case_when(
        difference_prevalence > 0 ~ "more prevalent in IBD",
        difference_prevalence < 0 ~ "less prevalent in IBD",
        TRUE ~ "no prevalence difference"
      )
    )

  ggplot(
    plot_df,
    aes(
      x = feature_short,
      y = difference_prevalence,
      fill = prevalence_direction
    )
  ) +
    geom_col() +
    coord_flip() +
    theme_minimal(base_size = 11) +
    labs(
      title = title,
      x = NULL,
      y = "Prevalence difference: IBD - healthy",
      fill = "Direction"
    )
}