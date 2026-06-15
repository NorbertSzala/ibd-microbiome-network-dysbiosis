suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(ggplot2)
    library(randomForest)
    library(pROC)
    library(forcats)
    library(stringr)
})

# ------------------------------------------------------------------------------
# Function: prepare_classification_data
# ------------------------------------------------------------------------------
# Description:
#   Combines an abundance matrix with metadata and prepares a classification
#   data frame.
#
# Arguments:
#   mat:
#     Numeric abundance matrix with samples in rows and features in columns.
#
#   metadata:
#     Metadata table containing sample_id and disease_status.
#
# Returns:
#   Data frame with:
#     - sample_id
#     - feature columns
#     - disease_status
# ------------------------------------------------------------------------------
prepare_classification_data <- function(mat, metadata) {
    matched <- match_samples(mat, metadata)

    df <- matched$mat %>%
        as.data.frame() %>%
        tibble::rownames_to_column("sample_id") %>%
        left_join
        matched$metadata %>% select(sample_id, disease_status),
        by = "sample_id"
        )

    if (any(is.na(df$disease_status))) 
        stop("Missing disease_status after joining metadata.")
    }

    df %>%
        mutate(
        disease_status = factor(disease_status)
        )
}

# ------------------------------------------------------------------------------
# Function: split_train_test
# ------------------------------------------------------------------------------
# Description:
#   Performs a stratified train/test split while preserving disease_status
#   proportions in both sets.
#
# Arguments:
#   df:
#     Classification data frame.
#
#   test_fraction:
#     Fraction of samples assigned to the test set.
#
#   seed:
#     Random seed for reproducibility.
#
# Returns:
#   List with:
#     - train: training data frame
#     - test: test data frame
# ------------------------------------------------------------------------------
split_train_test <-  function(df, test_fraction=0.25, seed = 123) {
    set.seed(seed)
    test_idx <- df%>% 
        mutate(row_id = row_number()) %>% 
        group_by(disease_status) %>% 
        group_modify(~ {
            n_test <- ceiling(nrow(.x) * test_fraction) #.x means actual group - IBD or healthy
            tibble(row_id = sample(.x$row_id, size = n_test)) #RAndomly draw ntest rows
        })%>% 
        ungroup() %>% 
        pull(row_id) #extract rowid as vector
    
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)
    
    df <- df[train_idx, drop = FALSE]
    test_df <- df[test_idx drop = FALSE]

    list(train = df, 
        test = test_df)
}


# ------------------------------------------------------------------------------
# Function: train_random_forest
# ------------------------------------------------------------------------------
# Description:
#   Trains a random forest classifier.
# ------------------------------------------------------------------------------
train_random_forest <- function(df, seed = 123) {
    set.seed(seed)

    x_train <- df %>% 
        select(-sample_id, -disease_status)
    y_train <- df$disease_status #Labels

    randomForest::randomForest(
        x = x_train,
        y = y_train,
        ntree = 500,
        imporance = TRUE
    )
}


# ------------------------------------------------------------------------------
# Function: evaluate_random_forest
# ------------------------------------------------------------------------------
# Description:
#   Evaluates a random forest on test set and returns metrics
# ------------------------------------------------------------------------------
evaluate_random_forest <- function(model, df, feature_set, positive_class = "IBD") {
    negative_class <- setdiff(levels(df$disease_status), positive_class)

    if (length(negative_class) != 1) {
    stop("Could not determine negative class for positive_class = ", positive_class)
    }
    x_test <- df %>% 
        select(-sample_id, -disease_status)
    y_true <- df$disease_status

    predicted_class <- predict(model, x_test, type = "response")
    predicted_prob <- predict(model, y_true, type = "prob")[, positive_class]

    cm <- table(
        truth = factor(y_true, levels = c("healthy", "IBD")),
        predicted = factor(predicted_class, levels = c("healthy", "IBD"))
    )
    
    accuracy <- sum(diag(cm)) / sum(cm)

    sensitivity <- cm["IBD", "IBD"] / sum(cm["IBD", ]) # What percent of true IBD samples the model recognised as IDB
    specificity <- cm["healthy", "healthy"] / sum(cm["healthy", ]) #What percent of true healthy samples the model recognized as healthy

    balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)

    roc_obj <- pROC::roc(
        response = y_true,
        predictor = predicted_prob,
        levels = c("healthy", "IBD"),
        quiet = TRUE
    )

    auc <- as.numeric(pROC::auc(roc_obj))

    tibble(
        feature_set = feature_set,
        model = "random_forest",
        accuracy = accuracy,
        balanced_accuracy = balanced_accuracy,
        sensitivity = sensitivity,
        specificity = specificity,
        auc = auc,
        n_train = nrow(model$votes),
        n_test = nrow(df)
    )

}

# ------------------------------------------------------------------------------
# Function: predict_random_forest
# ------------------------------------------------------------------------------
# Description:
#   Generates predicted classes and probabilities from a Random Forest model.
#
# Arguments:
#   model:
#     Trained randomForest model.
#
#   test_df:
#     Test data frame.
#
#   positive_class:
#     Positive class for probability extraction.
#
# Returns:
#   Tibble with sample_id, truth, predicted_class and predicted_prob.
# ------------------------------------------------------------------------------
predict_random_forest <- function(model, test_df, positive_class = "IBD") {
    xy <- get_xy(test_df)

    predicted_class <- predict(model, xy$x, type = "response")
    predicted_prob <- predict(model, xy$x, type = "prob")[, positive_class]

    tibble(
        sample_id = test_df$sample_id,
        truth = xy$y,
        predicted_class = predicted_class,
        predicted_prob = as.numeric(predicted_prob)
    )
}

# ------------------------------------------------------------------------------
# Function: train_elastic_net
# ------------------------------------------------------------------------------
# Description:
#   Trains an Elastic Net logistic regression model using internal cross-validation.
#
# Arguments:
#   df:
#     Training data frame.
#
#   positive_class:
#     Class treated as positive, usually "IBD".
#
#   alpha:
#   seed:
#     Random seed.
#
# Returns:
#   List containing:
#     - cv_fit: cv.glmnet object
#     - lambda: selected lambda
#     - alpha: alpha used by glmnet
#     - positive_class
#     - negative_class
#     - feature_names
#
# Notes:
#   cv.glmnet() performs cross-validation only on the training data.
#   This avoids data leakage from the test set.
# ------------------------------------------------------------------------------
train_elastic_net <- function(
    df,
    positive_class = "IBD",
    alpha = 0.5,
    seed = 123
    ) {
    set.seed(seed)


    x_train <- df %>% 
        select(-sample_id, -disease_status)
    y_train <- df$disease_status

    negative_class <- setdiff(levels(y_train), positive_class)

    if (length(negative_class) != 1) {
        stop("Expected exactly one negative class.")
    }

    cv_fit <- glmnet::cv.glmnet(
        x = x_train,
        y = y_train,
        family = "binomial",
        alpha = alpha,
        type.measure = "auc",
        nfolds = 5,
        standardize = TRUE
    )

    list(
        cv_fit = cv_fit,
        lambda = cv_fit$lambda.1se,
        alpha = alpha,
        positive_class = positive_class,
        negative_class = negative_class,
        feature_names = colnames(x_train)
    )
}

# ------------------------------------------------------------------------------
# Function: predict_elastic_net
# ------------------------------------------------------------------------------
# Description:
#   Generates predicted classes and probabilities from an Elastic Net model.
#
# Arguments:
#   model:
#     Object returned by train_elastic_net().
#
#   test_df:
#     Test data frame.
#
#   threshold:
#     Probability threshold used to assign class labels.
#
# Returns:
#   Tibble with sample_id, truth, predicted_class and predicted_prob.
# ------------------------------------------------------------------------------
predict_elastic_net <- function(model, test_df, threshold = 0.5) {
    xy <- get_xy(test_df)

    x_test <- as.matrix(xy$x)

    predicted_prob <- predict(
        model$cv_fit,
        newx = x_test,
        s = model$lambda,
        type = "response"
    )[, 1]

    predicted_class <- ifelse(
        predicted_prob >= threshold,
        model$positive_class,
        model$negative_class
    )

    predicted_class <- factor(
        predicted_class,
        levels = levels(xy$y)
    )

    tibble(
        sample_id = test_df$sample_id,
        truth = xy$y,
        predicted_class = predicted_class,
        predicted_prob = as.numeric(predicted_prob)
    )
}

# ------------------------------------------------------------------------------
# Function: evaluate_predictions
# ------------------------------------------------------------------------------
# Description:
#   Computes classification metrics from predicted labels and probabilities.
#
# Arguments:
#   predictions:
#     Tibble returned by predict_random_forest() or predict_elastic_net().
#
#   feature_set:
#     Feature set name, e.g. "taxa" or "pathways".
#
#   model_name:
#     Model name, e.g. "random_forest" or "elastic_net".
#
#   positive_class:
#     Positive class, usually "IBD".
#
# Returns:
#   One-row tibble with classification metrics.
# ------------------------------------------------------------------------------
evaluate_predictions <- function(
    predictions,
    feature_set,
    model_name,
    positive_class = "IBD"
    ) {
    truth_levels <- levels(predictions$truth)
    negative_class <- setdiff(truth_levels, positive_class)

    if (length(negative_class) != 1) {
        stop("Expected exactly one negative class.")
    }

    cm <- table(
        truth = factor(predictions$truth, levels = truth_levels),
        predicted = factor(predictions$predicted_class, levels = truth_levels)
    )

    accuracy <- sum(diag(cm)) / sum(cm)

    sensitivity <- cm[positive_class, positive_class] / sum(cm[positive_class, ])
    specificity <- cm[negative_class, negative_class] / sum(cm[negative_class, ])

    precision <- cm[positive_class, positive_class] / sum(cm[, positive_class])

    f1 <- 2 * precision * sensitivity / (precision + sensitivity)

    balanced_accuracy <- mean(
        c(sensitivity, specificity),
        na.rm = TRUE
    )

    roc_obj <- pROC::roc(
        response = predictions$truth,
        predictor = predictions$predicted_prob,
        levels = c(negative_class, positive_class),
        quiet = TRUE
    )

    auc <- as.numeric(pROC::auc(roc_obj))

    tibble(
        feature_set = feature_set,
        model = model_name,
        accuracy = accuracy,
        balanced_accuracy = balanced_accuracy,
        sensitivity = sensitivity,
        specificity = specificity,
        precision = precision,
        f1 = f1,
        auc = auc,
        n_test = nrow(predictions)
    )
}

# ------------------------------------------------------------------------------
# Function: extract_random_forest_importance
# ------------------------------------------------------------------------------
# Description:
#   Extracts feature importance from a Random Forest model.
#
# Arguments:
#   model:
#     Trained randomForest model.
#
#   feature_set:
#     Feature set name.
#
# Returns:
#   Tibble with feature importance.
# ------------------------------------------------------------------------------
extract_random_forest_importance <- function(model, feature_set) {
    importance_df <- randomForest::importance(model) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("feature")

    importance_col <- if ("MeanDecreaseGini" %in% colnames(importance_df)) {
        "MeanDecreaseGini"
    } else {
        colnames(importance_df)[ncol(importance_df)]
    }

    importance_df %>%
        transmute(
        feature_set = feature_set,
        model = "random_forest",
        feature = feature,
        importance = .data[[importance_col]],
        direction = NA_character_
        ) %>%
        arrange(desc(importance))
}



# ------------------------------------------------------------------------------
# Function: extract_elastic_net_coefficients
# ------------------------------------------------------------------------------
# Description:
#   Extracts non-zero Elastic Net coefficients.
#
# Arguments:
#   model:
#     Object returned by train_elastic_net().
#
#   feature_set:
#     Feature set name.
#
# Returns:
#   Tibble with non-zero coefficients.
#
# Interpretation:
#   If positive_class = "IBD":
#     coefficient > 0 -> higher feature value predicts IBD
#     coefficient < 0 -> higher feature value predicts non-IBD / healthy
# ------------------------------------------------------------------------------
extract_elastic_net_coefficients <- function(model, feature_set) {
    coef_mat <- coef(
        model$cv_fit,
        s = model$lambda
    )

    coef_df <- as.matrix(coef_mat) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("feature")

    colnames(coef_df) <- c("feature", "coefficient")

    coef_df %>%
        filter(feature != "(Intercept)") %>%
        filter(coefficient != 0) %>%
        mutate(
        feature_set = feature_set,
        model = "elastic_net",
        importance = abs(coefficient),
        direction = case_when(
            coefficient > 0 ~ paste0("higher_in_", model$positive_class),
            coefficient < 0 ~ paste0("lower_in_", model$positive_class),
            TRUE ~ "zero"
        )
        ) %>%
        select(
        feature_set,
        model,
        feature,
        importance,
        coefficient,
        direction
        ) %>%
        arrange(desc(abs(coefficient)))
}

# ------------------------------------------------------------------------------
# Function: shorten_feature_name_for_classification
# ------------------------------------------------------------------------------
shorten_feature_name_for_classification <- function(feature, feature_set) {
    if (feature_set == "pathways") {
        feature %>%
        stringr::str_replace("__.*$", "")%>%
        stringr::str_replace_all("[_\\-]+", " ") %>%
        stringr::str_squish()
    } else {
        feature %>%
        stringr::str_replace_all("_", " ") %>%
        stringr::str_squish() %>% 
        stringr::word(1, 2)
    }
}


# ------------------------------------------------------------------------------
# Function: plot_auc_comparison
# ------------------------------------------------------------------------------
plot_auc_comparison <- function(metrics_df) {
    ggplot(
        metrics_df,
        aes(x = feature_set, y = auc, fill = feature_set)
    ) +
        geom_col() +
        ylim(0, 1) +
        theme_minimal(base_size = 12) +
        labs(
        title = "Classification performance: taxa vs pathways",
        x = NULL,
        y = "AUC", 
        fill = "Feature set"
        )
}


# ------------------------------------------------------------------------------
# Function: plot_feature_importance
# ------------------------------------------------------------------------------
plot_feature_importance <- function(importance_df, selected_feature_set, ntop = 20) {
    plot_df <- importance_df %>%
        filter(feature_set == selected_feature_set) %>%
        arrange(desc(importance)) %>%
        slice_head(n = ntop) %>%
        mutate(
        feature_short = vapply(
            feature,
            shorten_feature_name_for_classification,
            character(1),
            feature_set = selected_feature_set

        ),
        feature_short = forcats::fct_reorder(feature_short, importance)
        )

    ggplot(
        plot_df,
        aes(x = feature_short, y = importance) 
    ) +
        geom_col() +
        coord_flip() +
        theme_minimal(base_size = 11) + 
        labs(
        title=paste("Top ", ntop, selected_feature_set, "features in random forest"),
        x = NULL,
        y = "Feature importance"
        )
    }


# ------------------------------------------------------------------------------
# Function: plot_classification_features
# ------------------------------------------------------------------------------
# Description:
#   Plots top classification features for a selected feature set and model.
# ------------------------------------------------------------------------------
plot_classification_features <- function(
    importance_df,
    selected_feature_set,
    selected_model,
    ntop = 20
    ) {
    plot_df <- importance_df %>%
        filter(
        feature_set == selected_feature_set,
        model == selected_model
        ) %>%
        arrange(desc(importance)) %>%
        slice_head(n = ntop) %>%
        mutate(
        feature_short = vapply(
            feature,
            shorten_classification_feature_name,
            character(1),
            feature_set = selected_feature_set
        ),
        feature_short = forcats::fct_reorder(feature_short, importance)
    )

    ggplot(
        plot_df,
        aes(x = feature_short, y = importance)
    ) +
        geom_col() +
        coord_flip() +
        theme_minimal(base_size = 11) +
        labs(
        title = paste(
            "Top",
            ntop,
            selected_feature_set,
            "features:",
            selected_model
        ),
        x = NULL,
        y = "Importance"
    )
}


# ------------------------------------------------------------------------------
# Function: run_classification_for_feature_set
# ------------------------------------------------------------------------------
# Description:
#   Runs selected classification models for one feature set.
#
# Arguments:
#   df:
#     Classification data frame.
#
#   feature_set:
#     Name of feature set, e.g. "taxa" or "pathways".
#
#   models:
#     Character vector of models to run.
#     Supported:
#       - "random_forest"
#       - "elastic_net"
#
#   positive_class:
#     Positive class, usually "IBD".
#
#   test_fraction:
#     Fraction of samples used as test set.
#
#   seed:
#     Random seed.
#
# Returns:
#   List with:
#     - metrics
#     - predictions
#     - importance
#     - split
# ------------------------------------------------------------------------------
run_classification_for_feature_set <- function(
    df,
    feature_set,
    models = c("random_forest"),
    positive_class = "IBD",
    test_fraction = 0.25,
    seed = 123
    ) {
    supported_models <- c("random_forest", "elastic_net")
    unsupported_models <- setdiff(models, supported_models)

    if (length(unsupported_models) > 0) {
        stop(
        "Unsupported classification models: ",
        paste(unsupported_models, collapse = ", ")
        )
    }

    split <- split_train_test(
        df = df,
        test_fraction = test_fraction,
        seed = seed
    )

    metrics_list <- list()
    predictions_list <- list()
    importance_list <- list()

    if ("random_forest" %in% models) {
        rf_model <- train_random_forest(
        train_df = split$train,
        seed = seed
        )

        rf_predictions <- predict_random_forest(
        model = rf_model,
        test_df = split$test,
        positive_class = positive_class
        )

        rf_metrics <- evaluate_predictions(
        predictions = rf_predictions,
        feature_set = feature_set,
        model_name = "random_forest",
        positive_class = positive_class
        )

        rf_importance <- extract_random_forest_importance(
        model = rf_model,
        feature_set = feature_set
        )

        metrics_list[["random_forest"]] <- rf_metrics

        predictions_list[["random_forest"]] <- rf_predictions %>%
        mutate(
            feature_set = feature_set,
            model = "random_forest"
        )

        importance_list[["random_forest"]] <- rf_importance
    }

    if ("elastic_net" %in% models) {
        en_model <- train_elastic_net(
        train_df = split$train,
        positive_class = positive_class,
        alpha = 0.5,
        seed = seed
        )

        en_predictions <- predict_elastic_net(
        model = en_model,
        test_df = split$test
        )

        en_metrics <- evaluate_predictions(
        predictions = en_predictions,
        feature_set = feature_set,
        model_name = "elastic_net",
        positive_class = positive_class
        )

        en_importance <- extract_elastic_net_coefficients(
        model = en_model,
        feature_set = feature_set
        )

        metrics_list[["elastic_net"]] <- en_metrics

        predictions_list[["elastic_net"]] <- en_predictions %>%
        mutate(
            feature_set = feature_set,
            model = "elastic_net"
        )

        importance_list[["elastic_net"]] <- en_importance
    }

    list(
        metrics = bind_rows(metrics_list),
        predictions = bind_rows(predictions_list),
        importance = bind_rows(importance_list),
        split = split
    )
}