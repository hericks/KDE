#' The S3 Class \code{Kernel} tries to ensure some of the properties of kernels
#' and is a subclass of the \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}}
#' (see: 'Details' for exact requirements).
#'
#' @details{
#' A kernel function is a real valued, integrable function, such that its integral over the real numbers equals one.
#' Kernel functions as \code{R} functions are required to:
#'
#'   1. be vectorised in its argument, taking a single numeric argument,
#'   returning a numerical vector of the same length only
#'
#'   2. return zero for inputs outside their support
#'
#'   3. can be integrated over their support using \code{integrate} without
#'   throwing an error
#'
#'   4. the integral over its support, using \code{integrate}, should evaluate to one.
#'
#'   The functions in this package don't just take \code{R} functions satisfying
#'   these conditions, but objects of S3 class \code{IntegrableFunction} (or one
#'   of its subclasses \code{Kernel}, \code{Density}).
#'
#'   The S3 class \code{Kernel} exists to ensure some of the most
#'   basic properties of kernel functions. The class is build on lists
#'   containing two named entries \code{fun} and \code{support}.
#'
#'   * \strong{`fun`} is an \code{R} function (the represented function) taking
#'   a single numeric argument and returning a numeric vector of the same
#'   length. This function should return zero outside of the interval given in
#'   the \code{support} entry. Using \code{integrate} over its support should evaluate to one.
#'
#'   * \strong{`support`} is a numeric vector of length 2 containing a lower-
#'   and upperbound for the support of the function stored in \code{fun} in its
#'   first and second entry respectively. In particular the values \code{-Inf}
#'   and \code{Inf} are allowed.
#'
#'   The constructor \code{Kernel} tries to construct a valid
#'   \code{Kernel} object based on the passed arguments. Returned
#'   objects are guaranteed to pass the validator [validate_Kernel]
#'   (\bold{Attention:} This does not guarantee the conditions in the first
#'   'Details' paragraph: see [validate_Kernel].).
#'
#' \describe{
#' \strong{List of built-in kernels functions:}
#'   \item{\code{\link[KDE:rectangular]{rectangular}}}
#'   \item{\code{\link[KDE:triangular]{triangular}}}
#'   \item{\code{\link[KDE:epanechnikov]{epanechnikov}}}
#'   \item{\code{\link[KDE:biweight]{biweight}}}
#'   \item{\code{\link[KDE:triweight]{triweight}}}
#'   \item{\code{\link[KDE:tricube]{tricube}}}
#'   \item{\code{\link[KDE:gaussian]{gaussian}}}
#'   \item{\code{\link[KDE:cosine]{cosine}}}
#'   \item{\code{\link[KDE:logistic]{logistic}}}
#'   \item{\code{\link[KDE:sigmoidFunction]{sigmoidFunction}}}
#'   \item{\code{\link[KDE:silverman]{silverman}}}
#'   }
#' }
#'
#' @param fun an \code{R} function taking a single numeric argument and
#' returning a numeric vector of the same length: see 'Details'.
#' @param support a numerical vector of length 2 containing the lower- and
#'   upperbound in the first and second entry respectively.
#'   \code{\link[KDE:Kernel]{Kernel}} will try to find bounds on the support itself if
#'   \code{NULL} is passed.
#'
#' @examples
#'
#'
#' @seealso
#' \code{\link[KDE:validate_Kernel]{validate_Kernel}}
#' \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}
#'
#' @export
Density <- function(fun, support = NULL){
  den <- new_Density(fun, support)
  validate_Density(den)
  den
}

#'
#'@include integrable_function.R
#'
#' @export
validate_Density <- function(x){
  stopifnot("object has to be of class Density"=inherits(x, "Density"))
  validate_IntegrableFunction(x)
  object <- x$fun
  lower <- x$support[1]
  upper <- x$support[2]

  testing_points <- seq(max(lower, -1e10), min(upper, 1e10), length.out=1e4)
  # prevent double checking with validate_integrable function
  stopifnot("density functions are non-negative" = all(object(testing_points) >= 0))

  stopifnot("The integral of a density over its support has to be one"=
              isTRUE(all.equal(integrate(object, lower = lower, upper = upper)[[1]], 1)))
  invisible(x)
}

new_Density <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Density"))
}
