#' Select Coefficients Based on Confidence Interval Bounds
#'
#' Given a data frame of bootstrap coefficients and a confidence level (1 - alpha),
#' this function returns the names of variables whose confidence intervals do not
#' include zero (i.e., either the lower bound > 0 or upper bound < 0).
#'
#' @param coef_df A data frame (or matrix) where each column contains bootstrap
#'     estimates for a single coefficient.
#' @param alpha A numeric value between 0 and 1 specifying the significance level
#'     (e.g., alpha = 0.05 for a 95% confidence interval).
#'
#' @return A character vector of selected coefficient names.
#' @export
select_predictors_by_ci <- function(coef_df, alpha) {
    ci_probs <- c(alpha / 2, 1 - alpha / 2)
    ci_bounds <- apply(coef_df, 2, quantile, probs = ci_probs)
    names(which(ci_bounds[1, ] > 0 | ci_bounds[2, ] < 0))
}
