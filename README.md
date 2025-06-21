
# glmnetstrap

<!-- badges: start -->
<!-- badges: end -->

**glmnetstrap** implements a bootstrap-based method for estimating predictor
significance in regularized regression models using `glmnet`, as described in
*Statistical Learning with Sparsity* by Hastie, Tibshirani, and Wainwright
([link](https://hastie.su.domains/StatLearnSparsity/)).

This approach allows you to assess variable stability across resampled datasets
and is particularly useful in high-dimensional settings with many correlated
predictors.


## Installation

You can install the development version of `glmnetstrap` from GitHub:

```r
remotes::install_github("cmatKhan/glmnetstrap")
```

All functions are fully documented. eg use `?bootstrap_cv_glmnet` after
installation.

## Example

```r
library(glmnetstrap)
library(readr)
library(dplyr)
library(stringr)
library(purrr)

# Load response and predictor data
response_files = list.files("~/htcf_ref/data/yeast_database_modelling/pull_data_20250325/response_frames",
                            full.names = TRUE)
names(response_files) = str_remove(basename(response_files), ".csv")

response_dfs = map(response_files, read_csv)
predictors = read_csv("~/htcf_ref/data/yeast_database_modelling/pull_data_20250325/predictors_normalized.csv")

# Construct input data for a specific regulator
create_input_df = function(regulator, response_dfs, predictors){
  tmp_response = response_dfs[[regulator]]
  colnames(tmp_response) = c('target_symbol', 'response')
  
  out_df = left_join(tmp_response, predictors) %>%
    filter(target_symbol != regulator)
  
  rownames(out_df) = out_df$target_symbol
  as.data.frame(out_df[, names(out_df) != "target_symbol"])
}

# Create interaction model formula
create_interaction_formula = function(ptf, mtf_list, response_var = "response", no_const_term = TRUE){
  interactors = paste(ptf, mtf_list, sep = ":")
  prefix = if (no_const_term) paste0(response_var, " ~ -1 + ", ptf, " + ")
           else paste0(response_var, " ~ ", ptf, " + ")
  as.formula(paste0(prefix, paste(interactors, collapse = " + ")))
}

# Regulator of interest
regulator = "ACA1"
x = create_input_df(regulator, response_dfs, predictors)
x_formula = create_interaction_formula(regulator, setdiff(colnames(x), c(regulator, "response")))
x_mm = model.matrix(x_formula, x)

# Define stratified CV fold classes (optional)
stratification_classification <- function(binding_vector,
                                          perturbation_vector,
                                          bins = c(0, 8, 64, 512, Inf),
                                          bin_by_binding_only = FALSE) {
  stopifnot(length(binding_vector) == length(perturbation_vector))
  stopifnot(is.numeric(binding_vector), is.numeric(perturbation_vector))
  if (length(bins) < 2) {
    warning("The number of bins is less than 2. Returning all ones.")
    return(rep(1L, length(binding_vector)))
  }

  # Create labels (1, 2, ..., n_bins)
  labels <- seq_len(length(bins) - 1)

  # Rank and bin binding
  binding_rank <- rank(-binding_vector, ties.method = "min")  # descending
  binding_bin <- cut(binding_rank, breaks = bins, labels = labels, right = TRUE)

  if (bin_by_binding_only) {
    return(as.integer(binding_bin))
  }

  # Rank and bin perturbation
  perturbation_rank <- rank(-perturbation_vector, ties.method = "min")  # descending
  perturbation_bin <- cut(perturbation_rank, breaks = bins, labels = labels, right = TRUE)

  # Combine into unique stratification classes
  binding_bin <- as.integer(binding_bin)
  perturbation_bin <- as.integer(perturbation_bin)
  
  # Transform the binding_bin by multiplying it by a scale factor
  # (the number of labels). This ensures that combining
  # (binding_bin - 1) * len(labels) + perturbation_bin
  # produces a unique label for each (binding_bin, perturbation_bin) pair.
  # For example, if binding_bin = 2 and perturbation_bin = 3, and there are
  # 4 labels, then the final label will be (2 - 1) * 4 + 3 = 7.
  # This allows for stratification across both variables simultaneously.
  return((binding_bin - 1L) * length(labels) + perturbation_bin)
}

fold_classes <- stratification_classification(
  binding_vector = x[[regulator]],
  perturbation_vector = x[["response"]],
  bins = c(0, 8, 64, 512, Inf),
  bin_by_binding_only = TRUE
)

# verify that the fold classes distribute observations as expected

# this is how folds are created in glmnetstrap (with k parameterized via
# nfolds)
foldid <- caret::createFolds(as.factor(fold_classes), k = 4, list = FALSE)

# view a table
table(fold_classes, foldid)

# or plot
x$foldid <- foldid
x$fold_class <- fold_classes

plt = ggplot(x, aes_string(x = regulator,
                           y = "response",
                           color = "factor(fold_class)")) +
  geom_point(alpha = 0.7, size = 1) +
  facet_wrap(~ foldid, ncol = 2) +
  labs(
    title = "Stratification class within CV folds",
    x = regulator,
    y = "Response",
    color = "Strat. class"
  ) +
  theme_minimal()


plotly::ggplotly(plt)

# Fit bootstrap models
res <- bootstrap_cv_glmnet(
  n_bootstraps = 100,
  x = x_mm,
  y = x[["response"]],
  fold_classes = fold_classes
)

# get coefficients which are significant at the 98% CI level
select_predictors_by_ci(extract_bootstrap_predictors(res), alpha = 0.02)
```

