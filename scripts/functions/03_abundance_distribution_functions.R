# ==============================================================================
# Helper functions
# ==============================================================================

# ------------------------------------------------------------------------------
# Function: matrix_to_long
# ------------------------------------------------------------------------------
# Description:
#   Converts wide sample by feature matrix into long format
#
# Arguments:
#   mat_df:
#       dataframe with 'sample_id' column and other feature columns
# Returns:
#   A long format df with columns: sample_id, feature, abundance, feature_set
# ------------------------------------------------------------------------------
matrix_to_long <- function(mat_df, feature_set) {
    mat_df %>% 
        pivot_longer(
            cols = -sample_id,
            names_to = "feature",
            values_to = "abundance"
        ) %>% 
        mutate(feature_set = feature_set)
}


# ------------------------------------------------------------------------------
# Function: create_abundance_density_plot
# ------------------------------------------------------------------------------
# Description:
#   Creates basic QC plots for processed taxonomic and pathway abundance matrices.
#   The input matrices are already filtered and transformed by the preprocessing step.
#   Therefore, this script plots the transformed abundance values directly.
#
# Arguments:
#   long_df:
#       Long format abundance table
#   title:
#       Plot title
# Returns:
#   ggplot object
# ------------------------------------------------------------------------------
create_abundance_density_plot <- function(long_df, title) {
    total_values <- sum(!is.na(long_df$abundance))
    n_zero = sum(long_df$abundance ==0, na.rm = TRUE)
    fraction_zero = n_zero/total_values

    zero_label = paste0(
        "Zero values removed: ",
        n_zero, " / ", total_values, "(", round(100*fraction_zero,1), "%)"
    )


    long_df %>% 
        filter(!is.na(abundance), abundance > 0) %>% 
        ggplot(aes(x=abundance)) + 
        geom_density(linewidth=0.8) +
        theme_minimal(base_size = 12) + 
        annotate("label", 
            x=Inf, y=Inf, label = zero_label, hjust = 1.05, vjust = 1.2, size = 3.5)+
        labs(title = title,
            x = 'Transformed abundance (log1p)',
            y="Density")
}

# ------------------------------------------------------------------------------
# Function: create_abundance_hist_plot
# ------------------------------------------------------------------------------
# Description:
#   Creates a histogram of log10-transformed non-zero abundance values. Zero values are removed before plotting, and the number and percentage of removed zero values are displayed as a plot annotation.
#
# Arguments:
#   long_df:
#       Long-format  table 

#   title:
#       Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
create_abundance_hist_plot <- function(long_df, title) {
    total_values <- sum(!is.na(long_df$abundance))
    n_zero <- sum(long_df$abundance == 0, na.rm = TRUE)
    fraction_zero <- n_zero / total_values

    zero_label <- paste0(
        "Zero values removed: ",
        n_zero, " / ", total_values, " (", round(100 * fraction_zero, 1),     "%)"
    )

    long_df %>%
        filter(!is.na(abundance), abundance > 0) %>%
        ggplot(aes(x = abundance)) +
        geom_histogram(aes(y=after_stat(count/sum(count))), bins = 100, alpha = 0.8) +
        annotate(
            "label", x = Inf, y = Inf,label = zero_label, hjust = 1.05, vjust = 1.2, size = 3.5
        ) +
        theme_minimal(base_size = 12) +
        labs(
            title = title,
            x = "Transformed abundance (log1p)",
            y = "Fraction of non-zero values"
        )
}


# ------------------------------------------------------------------------------
# Function: summarize_abundance
# ------------------------------------------------------------------------------
# Description:
#   Computes basic statistics for abudance values in given feature set
#
# Arguments:
#   long_df:
#       df abundance table
#
# Returns:
#   A table with abundance statistics
# ------------------------------------------------------------------------------
summarize_abundance <- function(long_df) {
    nonzero_abundance <- long_df$abundance[!is.na(long_df$abundance) & long_df$abundance > 0]

    if (length(nonzero_abundance) == 0) {
        nonzero_abundance <- NA_real_
    }

    long_df %>%
        summarize(
            feature_set = first(feature_set),

            n_values = sum(!is.na(abundance)),
            n_zero = sum(abundance == 0, na.rm = TRUE),
            n_nonzero = sum(abundance > 0, na.rm = TRUE),
            fraction_zero = n_zero / n_values,
            fraction_nonzero = n_nonzero / n_values,

            min_all = min(abundance, na.rm = TRUE),
            q25_all = quantile(abundance, 0.25, na.rm = TRUE),
            median_all = median(abundance, na.rm = TRUE),
            mean_all = mean(abundance, na.rm = TRUE),
            q75_all = quantile(abundance, 0.75, na.rm = TRUE),
            max_all = max(abundance, na.rm = TRUE),

            min_nonzero = min(nonzero_abundance, na.rm = TRUE),
            q25_nonzero = quantile(nonzero_abundance, 0.25, na.rm = TRUE),
            median_nonzero = median(nonzero_abundance, na.rm = TRUE),
            mean_nonzero = mean(nonzero_abundance, na.rm = TRUE),
            q75_nonzero = quantile(nonzero_abundance, 0.75, na.rm = TRUE),
            max_nonzero = max(nonzero_abundance, na.rm = TRUE)
        )
}
