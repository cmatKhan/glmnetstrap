#' Bootstrap Multiple cv.glmnet Models with Resampled Weights
#'
#' Performs repeated bootstrapped Lasso (or Elastic Net) modeling by generating
#' bootstrap weights and calling `fit_cv_glmnet()` for each replicate.
#'
#' @param n_bootstraps Number of bootstrap replicates.
#' @param seed Optional integer seed. When provided, ensures repeatability
#'     across replicates by offsetting the seed for each bootstrap
#'     (i.e., `seed + i`). If `NULL`, results are different each time.
#'
#' @inheritDotParams fit_cv_glmnet
#'
#' @return A list of `cv.glmnet` objects, one per bootstrap replicate.
#' @export
bootstrap_cv_glmnet = function(n_bootstraps, seed = NULL, ...){
    args = list(...)
    fits = list()

    # TODO: parallelize me!
    for (i in seq_len(n_bootstraps)){

        seed_i = if(!is.null(seed)) seed + i else NULL

        args$weights = generate_bootstrap_weights(length(args$y), seed = seed_i)

        fits[[i]] <- do.call(fit_cv_glmnet, args)
    }

    fits
}
