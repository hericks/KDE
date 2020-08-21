#' Kernels
#'
#' @description The kernels are key-objects of the \code{KDE} package. The S3
#'   class \code{Kernel} tries to ensure some of the properties of kernels and
#'   is a subclass of \code{\link{IntegrableFunction}} (see: 'Details' for exact
#'   requirements).
#'
#' @param fun a \code{R} function taking a single numeric argument and returning
#'   a numeric vector of the same length. See 'Details' for further
#'   requirements.
#' @param support numerical vector of length 2; the lower- and upper bound of
#'   the compact support in the first and second entry respectively. In
#'   particular non-finite values are prohibited. \code{IntegrableFunction} will
#'   try to find bounds on the support itself if \code{NULL} is passed.
#' @param subdivisions positive numeric scalar; the subdivisions parameter for
#'   the function \code{\link{integrate_primitive}}.
#' @param ... additional parameters to keep fixed during the evaluation of
#'   \code{fun}.
#'
#' @details A kernel function is a real valued, integrable function, such that
#'   its integral over the real numbers equals one. Kernel functions as \code{R}
#'   functions are required to
#'
#'   1. be vectorised in its argument, taking a single numeric argument,
#'   returning a numerical vector of the same length only,
#'
#'   2. return zero for inputs outside their compact support,
#'
#'   3. can be integrated over their support using \code{integrate_primitive}
#'   and the given number of subdivisions (the relative error converges).
#'
#'   4. yield an integral of nearly 1 (absolute error < 1%).
#'
#'   See the 'Details' section of \code{\link{IntegrableFunction}} for comments
#'   on the restrictiveness of compact supports.
#'
#'   The S3 class \code{Density} exists to ensure some of the most basic
#'   properties of density functions. The class is build on
#'   \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}} and inherits its
#'   structure.
#'
#'   The constructor \code{Kernel} tries to construct a valid \code{Kernel}
#'   object based on the passed arguments. Returned objects are guaranteed to
#'   pass the validator \code{\link{validate_Kernel}}.
#'
#'   \bold{Attention:} This does not guarantee the conditions in the first
#'   'Details' paragraph: see \code{validate_Kernel}.
#'
#'   \strong{List of built-in kernels functions:}
#'   \describe{
#'     \item{\code{\link{rectangular}}}{}
#'     \item{\code{\link{triangular}}}{}
#'     \item{\code{\link{epanechnikov}}}{}
#'     \item{\code{\link{biweight}}}{}
#'     \item{\code{\link{triweight}}}{}
#'     \item{\code{\link{tricube}}}{}
#'     \item{\code{\link{gaussian}}}{}
#'     \item{\code{\link{cosine}}}{}
#'     \item{\code{\link{logistic}}}{}
#'     \item{\code{\link{sigmoid}}}{}
#'     \item{\code{\link{silverman}}}{}
#'   }
#'
#'
#' @seealso \code{\link{validate_Kernel}} for the corresponding validator,
#' \code{\link{IntegrableFunction}} for more information about the superclass.
#'
#' @examples
#' rectangular_function <- function(u) 1/2*(abs(u) <= 1)
#' rectangular_ker <- Kernel(rectangular_function, support=c(-1,1))
#' x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, gaussian$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x, triangular$fun(x), col="red")
#' lines(x, rectangular_ker$fun(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "triangular","rectangular"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#'
#' @include integrable_function.R
#'
#' @export
Kernel <- function(fun, support = NULL, subdivisions=1000L, ...) {
  extra_args <- list(...)
  new_fun <- function(x) do.call(fun, c(list(x), extra_args))
  kern <- new_Kernel(new_fun, support, subdivisions)
  validate_Kernel(kern)
  kern
}

#' Validator for S3 class \code{Kernel}
#'
#' @description This function serves as a validator for the S3 class
#'   \code{Kernel}. See 'Details' for further information and
#'   potential flaws.
#'
#' @param x a \code{R} object to validate as an object of S3 class
#'   \code{\link{Kernel}}.
#'
#' @details The validator \code{validate_Kernel} can be used to verify objects
#'   as formally correct S3 objects of class \code{\link{Kernel}}. In particular
#'   the formal structure is ensured and
#'   \code{\link{validate_IntegrableFunction}} is called internally.
#'   Additionally this function \emph{tries to} (see 'Special Attention')
#'   validate the additional conditions of valid kernel functions (as specified
#'   in the first 'Details'-paragraph of \code{Kernel}).
#'
#' @section Special Attention:
#'
#'   Like all numerical routines, \code{validate_Kernel} can evaluate the
#'   represented function on a finite set of points. If the represented function
#'   returns valid results over nearly all its range, it is possible that the
#'   validator misses unexpected/wrong return values. Thus, using \code{Kernel}
#'   or \code{validate_Kernel} to construct and validate objects representing
#'   kernels is \emph{not} sufficient to ensure the properties listed in the
#'   first 'Details'-paragraph of \code{Kernel}, but serves more as a
#'   sanity-check.
#'
#' @seealso \code{\link{Kernel}} for more information about kernel functions and
#' the S3 class \code{Kernel}, \code{\link{IntegrableFunction}} for more
#' information about the superclass \code{IntegrableFunction}.
#'
#' @include integrable_function.R
#'
#' @export
validate_Kernel <- function(x){
  stopifnot("object has to be of class Kernel"=inherits(x, "Kernel"))
  validate_IntegrableFunction(x)
  object <- x$fun
  lower <- x$support[1]
  upper <- x$support[2]
  subdivisions = x$subdivisions
  stopifnot("The integral of a kernel over its support has to be one"=
              isTRUE(abs(integrate_primitive(object, lower = lower, upper = upper, subdivisions = subdivisions)$value-1) < 0.01))
  invisible(x)
}

new_Kernel <- function(fun, support, subdivisions=1000L, ..., subclass=NULL){
  new_IntegrableFunction(fun, support, subdivisions, subclass=c(subclass,"Kernel"))
}

#' Print objects of S3 class \code{Kernel}
#'
#' @param x object of S3 class \code{Kernel}; the object to print.
#' @param ... unused argument; used for compatibility with generic print.
#'
#' @export
print.Kernel <- function(x, ...) {
  print.IntegrableFunction(x, "Kernel")
}

