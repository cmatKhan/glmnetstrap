#' Fit a Cross-Validated glmnet Model with Optional Bootstrap Weights and
#'     Variable Subsetting
#'
#' This function fits a penalized regression model using `cv.glmnet` with
#'     support for: bootstrap-derived observation weights, stratified
#'     cross-validation folds, and an active set of coefficients. It also
#'     allows the user to override default settings such as `penalty.factor`
#'     and `weights` via `...`.
#'
#' @param x A numeric matrix of predictors (samples x features). Must have
#'     column names.
#' @param y A numeric response vector of length equal to the number of rows
#'     in `x`.
#' @param fold_classes Optional vector of class labels (same length as `y`)
#'     used to stratify the cross-validation folds. **NOTE**
#'     see \code{\link[caret]{createFolds}} with regards to how numeric vs
#'     factor levels are handled. If `NULL`, glmnet will handle CV splitting.
#'     If this is provided, `nfolds` must also be passed
#' @param active_predictors Optional logical or numeric index vector indicating
#'     which predictors to include in the model. If `NULL`, all predictors are
#'     used. If logical, must be the same length as `ncol(x)`. This argument
#'     is intended to be used in an interative process to drop out predictors
#' @param seed Integer seed for reproducible fold creation.
#'
#' @inheritDotParams glmnet::cv.glmnet
#'
#' @return A fitted `cv.glmnet` object.
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
#' y <- rnorm(100)
#' classes <- sample(c("A", "B"), 100, replace = TRUE)
#' fit <- cv_glmnet_stratified_active(x, y, fold_classes = classes, alpha = 1)
#' }
#'
#' @importFrom caret createFolds
#' @importFrom glmnet cv.glmnet
#' @export
fit_cv_glmnet <- function(
        x, y,
        fold_classes = NULL,
        active_predictors = NULL,
        seed = NULL,
        ...
) {

    # Extract optional arguments from ...
    args <- list(...)

    # If fold_classes provided, generate foldid
    if (!is.null(fold_classes) && is.null(args$nfolds)){
        stop("if `fold_classes` is provided, `nfolds` must also be specified")
    }
    else if (!is.null(fold_classes)) {
        if (!is.factor(fold_classes)){
            warning(paste0("`fold_classes` is not a factor. Check the ",
            "caret::createFolds() documentation to ensure that you ",
            "really want to use non-factored classes"))
        }
        set_optional_seed(seed)
        foldid <- caret::createFolds(as.factor(fold_classes),
                                     k = args$nfolds,
                                     list = FALSE)
        # Remove this from the args that will be passed into glmnet, since it
        # is handled by foldid
        args$nfolds <- NULL
    }

    if (is.null(active_predictors)) {
        active_predictors <- rep(TRUE, ncol(x))
    }

    # Ensure active_predictors is logical or numeric index
    if (is.logical(active_predictors) && length(active_predictors) != ncol(x)) {
        stop("Logical active_predictors must have length equal to ncol(x)")
    }

    x_internal <- x[, active_predictors, drop = FALSE]

    weights <- args$weights
    if (!is.null(weights) && length(weights) != length(y)) {
        stop("weights must be length equal to length(y)")
    }
    if (is.null(weights)) {
        weights <- rep(1, length(y))
    }

    penalty.factor <- args$penalty.factor
    if (is.null(penalty.factor)) {
        penalty.factor <- rep(1, ncol(x_internal))
    } else {
        if (length(penalty.factor) == ncol(x)) {
            penalty.factor <- penalty.factor[active_predictors]
        } else {
            stop("penalty.factor must be length equal to ncol(x)")
        }
    }

    # Remove handled arguments from args
    args$weights <- NULL
    args$penalty.factor <- NULL
    args$foldid <- NULL  # ignore user-supplied foldid

    # Build call
    base_args <- list(
        x = x_internal,
        y = y,
        weights = weights,
        penalty.factor = penalty.factor
    )

    if (!is.null(fold_classes)) {
        base_args$foldid <- foldid
    }

    do.call(glmnet::cv.glmnet, c(base_args, args))
}
