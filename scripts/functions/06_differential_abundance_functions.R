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

fisher_prevalence_test <- function(value, group) {
  df <- tibble(
    present = value > 0,
    group = group
  ) %>%
    filter(!is.na(present), !is.na(group)) %>%
    mutate(
      present = factor(present, levels = c(FALSE, TRUE)),
      group = factor(group, levels = c("healthy", "IBD"))
    )

  if (n_distinct(df$group) < 2) {
    return(NA_real_)
  }

  if (n_distinct(df$present) < 2) {
    return(NA_real_)
  }

  tab <- table(df$present, df$group)

  tryCatch(
    fisher.test(tab)$p.value,
    error = function(e) NA_real_
  )
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

      prevalence_healthy = mean(value[disease_status == "healthy"] > 0, na.rm = TRUE),
      prevalence_ibd = mean(value[disease_status == "IBD"] > 0, na.rm = TRUE),
      difference_prevalence = prevalence_ibd - prevalence_healthy,

      p_value = wilcoxon_test(
        value = value,
        group = disease_status
      ),
      p_value_prevalence = fisher_prevalence_test(
        value = value,
        group = disease_status),
      .groups = "drop"
    ) %>%
    group_by(feature_set) %>%
    mutate(
      p_adj = p.adjust(p_value, method = p_adjust_method),
      p_adj_prevalence = p.adjust(p_value_prevalence, method = p_adjust_method),
      direction = case_when(
        difference_mean > 0 ~ "higher_in_IBD",
        difference_mean < 0 ~ "lower_in_IBD",
        TRUE ~ "no_mean_difference"
      ),
      prevalence_direction = case_when(
      difference_prevalence > 0 ~ "more_prevalent_in_IBD",
      difference_prevalence < 0 ~ "less_prevalent_in_IBD",
      TRUE ~ "no_prevalence_difference"
    )) %>%
    ungroup() %>%
    arrange(feature_set, p_adj, desc(abs(difference_mean)))
}



# ------------------------------------------------------------------------------
# Function: get_top_features
# ------------------------------------------------------------------------------
# Description:
#  Selects top differential features based on adjusted p-value and absolute mean log1p abundance difference.
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
      desc(abs(difference_mean)),
      desc(abs(difference_prevalence))
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
  base_df <- results_df %>%
    filter(feature_set == selected_feature_set) %>%
    filter(!is.na(p_adj)) %>%
    filter(!is.na(difference_mean)) %>%
    filter(difference_mean != 0)

  plot_df <- bind_rows(
    base_df %>%
      filter(difference_mean > 0) %>%
      arrange(p_adj, desc(abs(difference_mean))) %>%
      slice_head(n = ceiling(ntop / 2)),

    base_df %>%
      filter(difference_mean < 0) %>%
      arrange(p_adj, desc(abs(difference_mean))) %>%
      slice_head(n = ceiling(ntop / 2))
  ) %>%
    mutate(
      feature_plot = forcats::fct_reorder(feature, difference_mean),
      direction = case_when(
        difference_mean > 0 ~ "higher in IBD",
        difference_mean < 0 ~ "lower in IBD",
        TRUE ~ "no mean difference"
      )
    )

  if (nrow(plot_df) == 0) {
    return(
      ggplot() +
        theme_void() +
        labs(title = paste(title, "- no abundance differences found"))
    )
  }

  ggplot(
    plot_df,
    aes(
      x = feature_plot,
      y = difference_mean,
      fill = direction
    )
  ) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(
      labels = function(x) {
        purrr::map_chr(
          x,
          ~ make_plot_label(.x, selected_feature_set)
        ) %>%
          stringr::str_wrap(
            width = ifelse(selected_feature_set == "pathways", 40, 28)
          )
      }
    ) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.001),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 8),

      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),

      axis.text.y = element_text(
        size = 7.8,
        lineheight = 0.82,
        margin = margin(r = 3)
      ),
      axis.text.x = element_text(size = 8.5),
      axis.title.x = element_text(size = 9.5),

      plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, margin = margin(b = 3)),
      plot.margin = margin(t = 2, r = 8, b = 2, l = 55)
    ) +
    labs(
      title = title,
      subtitle = "Top features in both directions; bars show mean log1p abundance difference",
      x = NULL,
      y = "Mean log1p abundance difference: IBD - healthy",
      fill = NULL
    )
}




plot_top_prevalence_features <- function(results_df, selected_feature_set, title, ntop = 20, max_label_words_taxa = 2, max_label_words_pathways = 3) {
  base_df <- results_df %>%
    filter(feature_set == selected_feature_set) %>%
    filter(!is.na(p_adj_prevalence)) %>%
    filter(!is.na(difference_prevalence)) %>%
    filter(difference_prevalence != 0)

  plot_df <- bind_rows(
    base_df %>%
      filter(difference_prevalence > 0) %>%
      arrange(p_adj_prevalence, desc(abs(difference_prevalence))) %>%
      slice_head(n = ceiling(ntop/2)),

    base_df %>%
      filter(difference_prevalence < 0) %>%
      arrange(p_adj_prevalence, desc(abs(difference_prevalence))) %>%
      slice_head(n = ceiling(ntop/2))
  ) %>%
    mutate(
      feature_plot = forcats::fct_reorder(feature, difference_prevalence),
      prevalence_direction = case_when(
        difference_prevalence > 0 ~ "more prevalent in IBD",
        difference_prevalence < 0 ~ "less prevalent in IBD",
        TRUE ~ "no prevalence difference"
      )
    )

  ggplot(
    plot_df,
    aes(
      x = feature_plot,
      y = difference_prevalence,
      fill = prevalence_direction
    )
  ) +
    geom_col(width = 0.72) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip(clip = "off") +
    scale_x_discrete(
      labels = function(x) {
        purrr::map_chr(
          x,
          ~ make_plot_label(.x, selected_feature_set)
        ) %>%
          stringr::str_wrap(
            width = ifelse(selected_feature_set == "pathways", 40, 28)
          )
      }
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0.04, 0.08))
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 8),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),

      axis.text.y = element_text(
        size = 7.8,
        lineheight = 0.82,
        margin = margin(r = 3)
      ),
      axis.text.x = element_text(size = 8.5),
      axis.title.x = element_text(size = 9.5),

      plot.title = element_text(size = 11, face = "bold", margin = margin(b = 2)),
      plot.subtitle = element_text(size = 9, margin = margin(b = 3)),
      plot.margin = margin(t = 2, r = 8, b = 2, l = 55)
    ) +
    labs(
      title = title,
      subtitle = "Features ranked by Fisher test on presence/absence",
      x = NULL,
      y = "Prevalence difference: IBD - healthy",
      fill = NULL
    )
}