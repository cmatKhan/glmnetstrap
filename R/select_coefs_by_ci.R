#' Select Predictors Whose Coefficients Are Consistently Nonzero
#'
#' This function selects predictors based on bootstrap coefficient
#'     distributions. For each predictor, it calculates a one-sided confidence
#'     interval (1 - \code{alpha}) and includes predictors whose bootstrap
#'     coefficient estimates consistently stay on one side of zero (i.e.,
#'     either the lower bound is strictly positive or the upper bound is
#'     strictly negative).
#'
#' This is equivalent to selecting predictors whose coefficients do not equal
#'     or cross zero more than \code{alpha} proportion of the time.
#'
#' @param coef_df A data frame or matrix where each column contains bootstrap
#'     estimates for a single coefficient.
#' @param alpha A numeric value between 0 and 1 specifying the exclusion
#'     threshold (e.g., \code{alpha = 0.05} excludes predictors whose
#'     coefficients cross zero in more than 5% of bootstrap replicates).
#'
#' @return A character vector of predictor names with stable, directionally
#'     consistent coefficients.
#' @export
select_predictors_by_ci <- function(coef_df, alpha) {
    ci_probs <- c(alpha / 2, 1 - alpha / 2)
    ci_bounds <- apply(coef_df, 2, quantile, probs = ci_probs)
    names(which(ci_bounds[1, ] > 0 | ci_bounds[2, ] < 0))
}
