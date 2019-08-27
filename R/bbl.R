#' @details
#' A typical workflow consists of the following steps:
#' \enumerate{
#'   \item Set up \code{bbl} object.
#'     Prepare a data frame containing data to be used as training set.
#'     Create a main object using data as input argument 
#'     (\code{\link{bbl}}).
#'   \item Train the model.
#'     See \code{\link{train}}.
#'     Perform cross-validation (\code{\link{crossval}}) 
#'     to optimize regularization.
#'   \item Make prediction on new data.
#'     See \code{\link{predict}}.
#' }
#' @keywords internal
#' @aliases bbl-package
"_PACKAGE"
#> [1] "_PACKAGE"