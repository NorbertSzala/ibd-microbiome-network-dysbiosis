suppressPackageStartupMessages({
    library(tidyverse)
    library(readr)
    library(ggplot2)
    library(vegan)
    library(stats)
})


# ------------------------------------------------------------------------------
# Function: remove_zero_variances
# ------------------------------------------------------------------------------
# Description:
#   Removes features with zero variances across samples (those features that does not change between samples). Those features doeas contribute to PCA
#
# Arguments:
#   mat:
#       numeric sample by feature matrix
#
# Returns:
#   Matrix without zer-variance columns
# ------------------------------------------------------------------------------
remove_zero_variances<- function(mat) {
    variances <- apply(mat, 2, var, na.rm=  TRUE) #Calculate variance for every column (feature). 
    keep <- variances > 0 & !is.na(variances)

    mat[, keep, drop = FALSE]
}


# ------------------------------------------------------------------------------
# Function: remove_empty_samples
# ------------------------------------------------------------------------------
# Description:
#   Removes samples with total abundance equals to zero. Those samples cant be used with bray curtis distance
# Arguments:
#   mat:
#     Numeric sample by feature matrix.
#
#   metadata:
#     Metadata table matched to matrix rows.
#
# Returns:
#   List with filtered matrix and metadata.
# ------------------------------------------------------------------------------
remove_empty_samples <- function(mat, metadata) {
  row_sums <- rowSums(mat, na.rm = TRUE)
  keep <- row_sums > 0 & !is.na(row_sums)

  n_removed <- sum(!keep)

  if (n_removed > 0) {
    message("Removing empty samples before Bray Curtis analysis: ", n_removed)
  }

  mat_filtered <- mat[keep, , drop = FALSE]
  metadata_filtered <- metadata[keep, , drop = FALSE]

  if (nrow(mat_filtered) == 0) {
    stop("All samples were removed as empty. Cannot continue ordination.")
  }

  list(
    mat = mat_filtered,
    metadata = metadata_filtered
  )
}

# ------------------------------------------------------------------------------
# Function: make_axis_limits
# ------------------------------------------------------------------------------
# Description:
#   Calculates axis limits using central quantiles. This removes the strongest
#   outliers from the displayed plotting range, but does not remove them from PCA,
#   PCoA, density contours, ellipses or centroid calculations.
#
# Arguments:
#   x:
#     Numeric vector with ordination coordinates.
#
#   probs:
#     Quantiles used to define the visible plotting range.
#
#   pad:
#     Small extra space added around the selected plotting range.
#
# Returns:
#   Numeric vector with lower and upper axis limit.
# ------------------------------------------------------------------------------
make_axis_limits <- function(x, probs = c(0.01, 0.99), pad = 0.02) {
  lim <- quantile(x, probs = probs, na.rm = TRUE)
  span <- diff(lim)

  if (!is.finite(span) || span == 0) {
    return(range(x, na.rm = TRUE))
  }

  c(lim[1] - pad * span, lim[2] + pad * span)
}

# ------------------------------------------------------------------------------
# Function: run_pca
# ------------------------------------------------------------------------------
# Description:
#   Runs PCA on the abundance matrix and returns sample coordinates and variance
#   explained by PCA axes. Features with zero variance and empty samples are
#   removed before PCA.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   metadata:
#     Metadata table matched to matrix rows. Must contain sample_id and
#     disease_status columns.
#
#   feature_set:
#     Name of feature set, for example "taxa" or "pathways".
#
# Returns:
#   List with PCA scores and variance explained table.
# ------------------------------------------------------------------------------
run_pca <- function(mat, metadata, feature_set, pca_scale = FALSE) {
  mat <- remove_zero_variances(mat)

  filtered <- remove_empty_samples(mat, metadata)
  mat <- filtered$mat
  metadata <- filtered$metadata

  if (any(is.na(mat))) {
    stop("Missing values detected in matrix before PCA.")
  }

  pca <- prcomp(mat, center = TRUE, scale. = pca_scale)

  explained <- (pca$sdev^2) / sum(pca$sdev^2)

  scores <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample_id") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(feature_set = feature_set)

  variance <- tibble(
    feature_set = feature_set,
    method = "PCA",
    axis = paste0("PC", seq_along(explained)),
    variance_explained = explained
  )

  list(
    scores = scores,
    variance = variance
  )
}

# ------------------------------------------------------------------------------
# Function: run_pcoa_bray
# ------------------------------------------------------------------------------
# Description:
#   Runs PCoA on Bray-Curtis distance matrix and returns sample coordinates and variance explained by ordination axes.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   metadata:
#     Metadata tibble containing `sample_id` and `disease_status`.
#
#   feature_set:
#       taxa/pathwyas for example
#
# Returns:
#  List with PCA scores and variance table

# I used Bray-Curtis as metric to count distances between each pair of samples. Then PCoA is made on this new matrix
# PCA = matrix (samples x feature)
# PCoA = matrix (samples x samples)
# ------------------------------------------------------------------------------
run_pcoa_bray <- function(mat, metadata, feature_set) {
  mat <- remove_zero_variances(mat)

  filtered <- remove_empty_samples(mat, metadata)
  mat <- filtered$mat
  metadata <- filtered$metadata

  if (any(mat < 0, na.rm = TRUE)) {
    stop(
      "Negative values detected. Bray-Curtis distance requires non-negative data. ",
      "Use log1p-transformed data, not CLR-transformed data, for Bray-Curtis PCoA."
    )
  }

  bray_dist <- vegan::vegdist(mat, method = "bray")

  if (any(is.na(bray_dist))) {
    stop(
      "Bray-Curtis distance contains NA values after removing empty samples. ",
      "Check whether the matrix contains missing values or invalid rows."
    )
  }

  pcoa <- cmdscale(
    bray_dist,
    eig = TRUE,
    k = 2
  )

  scores <- as.data.frame(pcoa$points) %>%
    setNames(c("PCoA1", "PCoA2")) %>%
    rownames_to_column("sample_id") %>%
    left_join(metadata, by = "sample_id") %>%
    mutate(feature_set = feature_set)

  positive_eig <- pcoa$eig[pcoa$eig > 0]
  variance_explained <- positive_eig / sum(positive_eig)

  variance <- tibble(
    feature_set = feature_set,
    method = "PCoA_BrayCurtis",
    axis = paste0("PCoA", seq_along(variance_explained)),
    variance_explained = variance_explained
  )

  list(
    scores = scores,
    variance = variance
  )
}

# ------------------------------------------------------------------------------
# Function: plot_ordination_clean
# ------------------------------------------------------------------------------
# Description:
#   Creates PCA/PCoA plot for presentation. Points are downsampled only
#   for readability, while density contours, ellipses and centroids are calculated
#   from all samples.
#
# Arguments:
#   scores:
#     Data frame with ordination scores and disease_status column.
#
#   x_col, y_col:
#     Names of columns used as x and y axis.
#
#   x_label, y_label:
#     Axis labels with explained variance.
#
#   title:
#     Plot title.
#
#   n_points_per_group:
#     Number of background points sampled from each disease group.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_ordination_clean <- function(
  scores,
  x_col,
  y_col,
  x_label,
  y_label,
  title,
  n_points_per_group = 500
) {
  group_colors <- c(
    "healthy" = "#4C72B0",
    "IBD" = "#C44E52"
  )

  plot_df <- scores %>%
    filter(
      is.finite(.data[[x_col]]),
      is.finite(.data[[y_col]]),
      !is.na(disease_status)
    ) %>%
    mutate(
      disease_status = factor(disease_status, levels = c("healthy", "IBD"))
    )

  # Points are downsampled only for plotting.
  points_df <- plot_df %>%
    group_by(disease_status) %>%
    group_modify(~ slice_sample(.x, n = min(nrow(.x), n_points_per_group))) %>%
    ungroup()

  # Centroids are calculated from all samples in each group.
  centroids <- plot_df %>%
    group_by(disease_status) %>%
    summarize(
      x_centroid = mean(.data[[x_col]], na.rm = TRUE),
      y_centroid = mean(.data[[y_col]], na.rm = TRUE),
      .groups = "drop"
    )

  # Axis limits are based on central 98% of points.
  x_lim <- make_axis_limits(plot_df[[x_col]], probs = c(0.01, 0.99), pad = 0.02)
  y_lim <- make_axis_limits(plot_df[[y_col]], probs = c(0.01, 0.99), pad = 0.02)

  ggplot(
    plot_df,
    aes(x = .data[[x_col]], y = .data[[y_col]], color = disease_status)
  ) +
    geom_point(
      data = points_df,
      alpha = 0.3,
      size = 0.75,
      stroke = 0
    ) +
    stat_density_2d(
      aes(group = disease_status),
      bins = 5,
      linewidth = 0.5,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    stat_ellipse(
      aes(group = disease_status),
      type = "norm",
      level = 0.7,
      linewidth = 0.75,
      show.legend = FALSE
    ) +
    geom_point(
      data = centroids,
      aes(
        x = x_centroid,
        y = y_centroid,
        color = disease_status
      ),
      inherit.aes = FALSE,
      shape = 4,
      size = 3.5,
      stroke = 1.2,
      show.legend = FALSE
    ) +
    coord_cartesian(
      xlim = x_lim,
      ylim = y_lim,
      expand = FALSE
    ) +
    scale_color_manual(values = group_colors) +
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          size = 2.5,
          alpha = 1
        )
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),
      legend.box = "horizontal",

      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),

      axis.title = element_text(size = 10),
      axis.text = element_text(size = 8),

      panel.grid.minor = element_blank(),

      plot.margin = ggplot2::margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1,
        unit = "cm"
      )
    ) +
    labs(
      title = title,
      x = x_label,
      y = y_label,
      color = "Disease status"
    )
}




# ------------------------------------------------------------------------------
# Function: plot_pca
# ------------------------------------------------------------------------------
# Description:
#   Creates PCA plot using cleaned ordination visualization.
#
# Arguments:
#   pca_result:
#     Output list returned by run_pca().
#
#   title:
#     Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_pca <- function(pca_result, title) {
  var_pc1 <- pca_result$variance %>%
    filter(axis == "PC1") %>%
    pull(variance_explained) %>%
    first()

  var_pc2 <- pca_result$variance %>%
    filter(axis == "PC2") %>%
    pull(variance_explained) %>%
    first()

  plot_ordination_clean(
    scores = pca_result$scores,
    x_col = "PC1",
    y_col = "PC2",
    x_label = paste0("PC1 (", round(var_pc1 * 100, 1), "%)"),
    y_label = paste0("PC2 (", round(var_pc2 * 100, 1), "%)"),
    title = title
  )
}


# ------------------------------------------------------------------------------
# Function: plot_pcoa
# ------------------------------------------------------------------------------
# Description:
#   Creates PCoA plot using cleaned ordination visualization.
#
# Arguments:
#   pcoa_result:
#     Output list returned by run_pcoa_bray().
#
#   title:
#     Plot title.
#
# Returns:
#   ggplot object.
# ------------------------------------------------------------------------------
plot_pcoa <- function(pcoa_result, title) {
  var_axis1 <- pcoa_result$variance %>%
    filter(axis == "PCoA1") %>%
    pull(variance_explained) %>%
    first()

  var_axis2 <- pcoa_result$variance %>%
    filter(axis == "PCoA2") %>%
    pull(variance_explained) %>%
    first()

  plot_ordination_clean(
    scores = pcoa_result$scores,
    x_col = "PCoA1",
    y_col = "PCoA2",
    x_label = paste0("PCoA1 (", round(var_axis1 * 100, 1), "%)"),
    y_label = paste0("PCoA2 (", round(var_axis2 * 100, 1), "%)"),
    title = title
  )
}

# ------------------------------------------------------------------------------
# Function: run_permanova
# ------------------------------------------------------------------------------
# Description:
#   Runs permanova using Bray-Curtis distance to test whether disease status explains global differences in microbiome composition.
#
# Arguments:
#   mat:
#     Numeric sample-by-feature matrix.
#
#   metadata:
#     Metadata tibble containing `disease_status`.
#
#   feature_set:
#     taxa/pathways
#
# Returns:
#   Tidy tibble with permanova results.
# ------------------------------------------------------------------------------
run_permanova <- function(
  mat,
  metadata,
  feature_set,
  group_col = "disease_status",
  permutations = 999,
  distance_method = "bray",
  seed = NULL
) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  if (!group_col %in% colnames(metadata)) {
    stop("Grouping column not found in metadata: ", group_col)
  }

  if (is.null(rownames(mat))) {
    stop("Matrix must have sample IDs as rownames.")
  }

  if (!"sample_id" %in% colnames(metadata)) {
    stop("Metadata must contain `sample_id` column.")
  }

  # Reorder metadata to matrix order.
  metadata <- metadata %>%
    filter(.data$sample_id %in% rownames(mat)) %>%
    arrange(match(.data$sample_id, rownames(mat)))

  mat <- mat[metadata$sample_id, , drop = FALSE]

  if (!identical(rownames(mat), metadata$sample_id)) {
    stop("Matrix rows and metadata sample_id are not aligned.")
  }

  # Keep only samples with non-missing group labels.
  keep_group <- !is.na(metadata[[group_col]])
  mat <- mat[keep_group, , drop = FALSE]
  metadata <- metadata[keep_group, , drop = FALSE]

  metadata[[group_col]] <- factor(metadata[[group_col]])

  group_counts <- table(metadata[[group_col]])

  if (length(group_counts) < 2) {
    stop("PERMANOVA requires at least two groups.")
  }

  if (any(group_counts < 2)) {
    stop(
      "Each group should contain at least 2 samples. Group counts: ",
      paste(names(group_counts), group_counts, sep = "=", collapse = ", ")
    )
  }

  # Remove features that cannot contribute to distances or PCA-like structure.
  mat <- remove_zero_variances(mat)

  if (ncol(mat) == 0) {
    stop("No features left after removing zero-variance features.")
  }

  # Remove samples with zero total abundance.
  filtered <- remove_empty_samples(mat, metadata)
  mat <- filtered$mat
  metadata <- filtered$metadata

  if (anyNA(mat)) {
    stop("Missing values detected in matrix before PERMANOVA.")
  }

  if (any(mat < 0, na.rm = TRUE)) {
    stop(
      "Negative values detected. Bray-Curtis PERMANOVA requires non-negative data. ",
      "Use log1p-transformed or untransformed non-negative abundance data, not CLR."
    )
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  dist_mat <- vegan::vegdist(mat, method = distance_method)

  if (any(is.na(dist_mat))) {
    stop(
      "Distance matrix contains NA values. Check empty samples, missing values, ",
      "or invalid abundance values."
    )
  }

  # PERMANOVA: tests whether group explains global distance structure.
  formula_text <- paste("dist_mat ~", group_col)
  permanova <- vegan::adonis2(
    as.formula(formula_text),
    data = metadata,
    permutations = permutations,
    by = "terms"
  )

  permanova_tbl <- as.data.frame(permanova) %>%
    rownames_to_column("term") %>%
    filter(term == group_col) %>%
    mutate(
      feature_set = feature_set,
      distance = distance_method,
      test = "PERMANOVA",
      n_samples = nrow(metadata),
      n_features = ncol(mat),
      permutations = permutations
    ) %>%
    select(
      feature_set,
      test,
      distance,
      term,
      n_samples,
      n_features,
      permutations,
      Df,
      SumOfSqs,
      R2,
      F,
      `Pr(>F)`
    ) %>%
    rename(p_value = `Pr(>F)`)

  # Dispersion test: checks whether groups differ in within-group variability.
  dispersion <- vegan::betadisper(
    dist_mat,
    group = metadata[[group_col]]
  )

  dispersion_test <- vegan::permutest(
    dispersion,
    permutations = permutations
  )

  dispersion_tbl <- tibble(
    feature_set = feature_set,
    test = "PERMDISP",
    distance = distance_method,
    term = group_col,
    n_samples = nrow(metadata),
    n_features = ncol(mat),
    permutations = permutations,
    Df = NA_real_,
    SumOfSqs = NA_real_,
    R2 = NA_real_,
    F = as.numeric(dispersion_test$tab[1, "F"]),
    p_value = as.numeric(dispersion_test$tab[1, "Pr(>F)"])
  )

  bind_rows(permanova_tbl, dispersion_tbl)
}