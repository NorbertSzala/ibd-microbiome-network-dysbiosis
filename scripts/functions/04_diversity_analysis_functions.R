# ==============================================================================
# 3. Helper functions
# ==============================================================================

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
calculate_diversity_metrics_vegan<- function(mat, feature_set, calculate_chao1=params$calculate_chao1) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  if (any(mat<0, na.rm = TRUE)) {
    stop(
      "Negative values detected in abundance matrix. Diversity metrics should be calculated only on row filtered non negative data"
    )
  }

  diversity_wide <- tibble(
    sample_id = rownames(mat),
    richness  = vegan::specnumber(mat), # typical richness, but counts even rarely 'organisms' or artefacts - detected taxa/paths
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
    filter(is.finite(value)) %>%
    summarize(
      n_healthy = sum(disease_status =='healthy'),
      n_IBD = sum(disease_status == 'IBD'),
      median_healthy = median(value[disease_status == 'healthy'], na.rm = TRUE),
      median_IBD = median(value[disease_status == 'IBD'], na.rm = TRUE),
      difference_median = median_IBD - median_healthy,
      p_value = if (n_healthy >= 2 && n_IBD >= 2) {
        wilcox.test(
          value ~ disease_status, #values are defined earlier metrics. It cause that the test will be made on each pairs (taxa/pathways + metrics),
            exact = FALSE
        )$p.value
      } else {
          NA_real_
      },

        .groups = "drop"
      ) %>% mutate (
        p_adj = p.adjust(p_value, method = "BH") #reduce false discovery rate
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