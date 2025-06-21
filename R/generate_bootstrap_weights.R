#' Generate Bootstrap Weights via Sampling
#'
#' Creates bootstrap weights for observations by sampling indices with replacement.
#' The number of times each observation is sampled becomes its weight.
#'
#' @param n Number of observations in the dataset.
#' @param seed Optional integer seed. Controls reproducibility of the resampling.
#'             If `NULL`, results will vary on each call.
#' @param standardize Logical; if `TRUE`, scales the weights so they sum to 1.
#'
#' @return Numeric vector of length `n` with bootstrap weights (integer if not standardized).
#' @export
generate_bootstrap_weights <- function(n, seed = NULL, standardize = FALSE) {
    set_optional_seed(seed)

    # Sample indices with replacement
    idx <- sample(seq_len(n), size = n, replace = TRUE)

    # Tabulate frequency of each observation in the resample
    w <- tabulate(idx, nbins = n)

    if (standardize) {
        w <- w / sum(w)
    }

    return(w)
}
