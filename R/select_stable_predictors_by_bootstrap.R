#' Iterative Bootstrap CV with Confidence Interval-Based Variable Selection
#'
#' This function performs iterative bootstrapped Lasso modeling using
#' \code{cv.glmnet()} across increasing confidence intervals (CIs),
#' dropping variables whose confidence intervals include zero. Iteration
#' stops when the number of selected variables stabilizes between CI steps.
#' A final bootstrap run is performed at the stabilized variable set.
#'
#' @param x A numeric matrix of predictors (samples × features).
#' @param y A numeric response vector of length equal to \code{nrow(x)}.
#' @param fold_classes A factor or vector used to create stratified folds.
#' @param active_predictors Optional logical or numeric index indicating
#'     which predictors to include. Default is all.
#' @param ci_start Starting confidence interval for iterative selection
#'     (default: 50).
#' @param ci_final Final confidence interval used in final modeling step
#'     (default: 98).
#' @param ci_step Step size for increasing CI between iterations (default: 2).
#' @param n_bootstraps Number of bootstrap replicates to run per CI level.
#' @param seed Random seed for reproducibility. Used for bootstrap sampling
#'     and fold generation.
#'
#' @inheritDotParams bootstrap_cv_glmnet
#'
#' @return A list containing:
#' \describe{
#'   \item{iterative_final_ci}{The final confidence interval at which
#'       stabilization occurred.}
#'   \item{stabilized_variables}{Character vector of variables selected
#'       after stabilization.}
#'   \item{fits}{A list of fitted \code{cv.glmnet} objects from the final
#'       bootstrap pass.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' x <- matrix(rnorm(100 * 10), ncol = 10)
#' y <- rnorm(100)
#' folds <- sample(rep(1:4, length.out = 100))
#' result <- select_stable_predictors_by_bootstrap(x, y, folds, alpha = 1)
#' result$stabilized_variables
#' }
#'
#' @importFrom futile.logger flog.debug
#' @export
select_stable_predictors_by_bootstrap <- function(
        x, y,
        fold_classes,
        ci_start = 50,
        ci_final = 98,
        ci_step = 2,
        n_bootstraps = 100,
        seed = NULL,
        ...
) {
    current_ci <- ci_start
    prev_n_vars <- NA_integer_
    stabilized_vars <- NULL
    iter <- 0

    # there are two conditions which break this loop -- the first is if there
    # are no selected coefficients. The second is if the number of selected model
    # coefficients is equal to the previous iteration.
    repeat {

        # execute the bootstrapping on the lasso model. Note that the `...`
        # args are passed through to glmnet.
        fits <- bootstrap_cv_glmnet(
            n_bootstraps = n_bootstraps,
            x = x,
            y = y,
            fold_classes = fold_classes,
            active_predictors = active_predictors,
            seed = seed,
            ...
        )

        coef_df <- extract_bootstrap_predictors(fits)

        selected = select_predictors_by_ci(coef_df, 1 - current_ci / 100)

        if (length(selected) == 0) {
            warning("No variables survive the iterative CI dropout")
            return(NULL)
        }

        futile.logger::flog.debug(
            "CI=%.1f selected variables [%s]",
            current_ci,
            paste(shQuote(selected), collapse = ", ")
        )

        n_selected <- length(selected)

        if (!is.na(prev_n_vars) && n_selected == prev_n_vars) {
            stabilized_vars <- selected
            break
        }

        rm(fits)
        gc()
        active_predictors <- colnames(x) %in% selected
        prev_n_vars <- n_selected
        current_ci <- min(current_ci + ci_step, ci_final)
        iter <- iter + 1
    }

    # Final pass with stabilized variable set and full CI
    final_fits <- bootstrap_cv_glmnet(
        n_bootstraps = n_bootstraps,
        x = x,
        y = y,
        fold_classes = fold_classes,
        seed = if(!is.null(seed)) seed + 100 else NULL,
        ...
    )

    final_coef_df <- extract_bootstrap_predictors(final_fits)
    final_selection_freqs <- selection_frequency(final_coef_df)

    list(
        iterative_final_ci = current_ci,
        input_coef_subset = stabilized_vars,
        fits = final_fits
    )
}
