#' Kernel
#'
#' @description
#' The kernels are the core of the KDE package.
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
#' rectangular_function <- function(u){
#'   check_kernel_conditions(u)
#'   return(1/2*(abs(u) <= 1))
#' }
#' rectangular_ker <- Kernel(rectangular_function, support=c(-1,1))
#' x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, gaussian$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x,triangular$fun(x), col="red")
#' lines(x, rectangular_ker$fun(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "triangular","rectangular"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#' @seealso
#' \code{\link[KDE:validate_Kernel]{validate_Kernel}}
#' \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}
#'
#' @export
Kernel <- function(fun, support = NULL) {
  kern <- new_Kernel(fun, support)
  validate_Kernel(kern)
  kern
}

#' Validator for S3 class \code{Kernel}
#'
#' @description This function serves as a validator for the S3 class
#'   \code{Kernel}. See 'Details' for further information and
#'   potential flaws.
#'
#' @param x an \code{R} object to validate as object of S3 class
#'   \code{[Kernel]}.
#'
#' @details The validator \code{validate_Kernel} can be used to
#'   verify objects as formally correct S3 objects of class
#'   [Kernel]. In particular the formal structure is ensured and it makes use of the [validate_IntegrableFunction.
#'   Additionally this function \emph{tries to} (see 'Special Attention')
#'   validate the additional conditions of valid integrable functions (as
#'   specified in the first 'Details'-paragraph of [Kernel]).
#'
#' @section Special Attention:
#'
#'   Like all numerical routines, \code{validate_Kernel} can
#'   evaluate the represented function on a finite set of points only. If the
#'   represented function returns valid results over nearly all its range, it is
#'   possible that this function misses unexpected/wrong return values. Thus,
#'   using [Kernel] or [validate_Kernel] to construct
#'   and validate objects representing integrable functions is \emph{not}
#'   sufficient to ensure the properties \[1-4\] listed in the first
#'   'Details'-paragraph of [Kernel], but serves more as a
#'   sanity-check.
#'
#' @seealso
#' * [Kernel] for more information about kernel functions and the S3 class \code{Kernel}.
#' * [IntegrableFunction] for more information about integrable functions and the S3 class \code{IntegrableFunctions}.
#' * [validate_IntegrableFunction] for more information about the IntegrableFunction validator.
#'
#'@export
validate_Kernel <- function(x){
  stopifnot("object has to be of class Kernel"=inherits(x, "Kernel"))
  validate_IntegrableFunction(x)
  object <- x$fun
  lower <- x$support[1]
  upper <- x$support[2]
  stopifnot("The integral of a density over its support has to be one"=
              isTRUE(all.equal(integrate(object, lower = lower, upper = upper)[[1]], 1)))
  invisible(x)
}

#' @include integrable_function.R
new_Kernel <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}
