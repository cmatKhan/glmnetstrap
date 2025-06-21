#' Optionally Set Random Seed
#'
#' Sets the random seed if a seed is provided. This function enables conditional
#' reproducibility—when `seed` is `NULL`, no seed is set and results may vary
#' between runs; otherwise, a fixed seed ensures reproducible outcomes.
#'
#' @param seed Optional integer seed. If `NULL`, randomness is not fixed.
#'
#' @return None. This function is called for its side effect.
#' @keywords internal
#' @export
set_optional_seed <- function(seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
}
