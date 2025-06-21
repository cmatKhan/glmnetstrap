#' Extract Coefficients from Bootstrap-Replicated Models
#'
#' Pulls the non-intercept coefficients at a specified regularization strength
#' from each `cv.glmnet` object in the input list.
#'
#' @param fits A list of fitted `cv.glmnet` objects.
#' @param s The lambda value at which to extract coefficients (default: `"lambda.min"`).
#'
#' @return A data frame of coefficients (rows = bootstraps, columns = features).
#' @export
extract_bootstrap_predictors <- function(fits, s = "lambda.min") {
    coef_results = lapply(fits, function(fit) {
        predictors <- as.vector(coef(fit, s = s))[-1]
        names(predictors) <- rownames(coef(fit))[-1]
        predictors
    })

    as.data.frame(do.call(rbind, coef_results))

}

#' Extract Intercepts from Bootstrap-Replicated Models
#'
#' @param fits A list of fitted `cv.glmnet` objects.
#'
#' @return A numeric vector of intercepts (one per bootstrap).
#' @export
extract_intercepts <- function(fits) {
    sapply(fits, function(fit) coef(fit, s = "lambda.min")[1])
}

#' Extract Optimal Lambda Values from Bootstrap Models
#'
#' @param fits A list of fitted `cv.glmnet` objects.
#'
#' @return A numeric vector of optimal lambda values (i.e., `lambda.min` from each fit).
#' @export
extract_lambda_min <- function(fits) {
    sapply(fits, function(fit) fit$lambda.min)
}

#' Compute Variable Selection Frequencies Across Bootstraps
#'
#' Calculates how often each coefficient is non-zero across a list or matrix of
#' coefficient vectors.
#'
#' @param fits A list of fitted `cv.glmnet` objects.
#'
#' @return A named numeric vector with selection frequencies (values between 0 and 1).
#' @export
selection_frequency <- function(fits) {
    colMeans(extract_bootstrap_predictors(fits) != 0)
}
