#' Kernel
#'
#' @description
#' The kernels in this package are S3 objects based on a list, holding a function and its support.
#' They are a subclass of the \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}}.
#'
#' @details{
#' A kernel is a real valued, integrable function, such that its integral over the real numbers equals one.
#' The built-in kernel functions are vectorised.
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
#'   returning a numeric vector of the same length: see 'Details'.
#' @param support a numerical vector of length 2 containing the lower- and
#'   upperbound in the first and second entry respectively.
#'   \code{IntegrableFunction} will try to find bounds on the support itself if
#'   \code{NULL} is passed: see 'Details' of \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}.
#'
#' TODO: more examples (integrate/support ins Spiel bringen)
#' @examples
#' x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, gaussian$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x,triangular$fun(x), col="red")
#' lines(x, rectangular$fun(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "triangular","rectangular"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#' @seealso
#' \code{\link[KDE:validate_Kernel]{validate_Kernel}}
#' \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}
#'
#' @export
Kernel <- function(fun, support) {
  kern <- new_Kernel(fun, support)
  validate_Kernel(kern)
  kern
}

#' Validation if object is a kernel function.
#' TODO: umschreiben/mit Kernel zusammenfassen
#'
#' @description
#' The \code{is_kernel()} function is used for validating a \link[KDE:Kernel]{Kernel}.
#' It has to be a integrable, numerical function with integral over R equal to one.
#'
#' @details
#' Hence a Kernel is a subclass of the \code{\link[KDE:IntegrableFunction]{IntegrableFunctions}},
#' it makes use of [validate_IntegrableFunction].
#' In addition, it will check wether the object is of class \link[KDE:Kernel]{Kernel} and if the integral over its support equals one.
#'
#' @param object the object to be tested.
#'
#' @return Boolean. Returns True if object is a kernel function, FALSE if not. Is not supossed to return an error.
#'
#' @section Issue:
#' Because of the use of the Base-R function integrate, functions that are almost everywhere zero can be kernels, but will not be detected.
#' For example: kernels that are shifted with a very small bandwidth.
#'
#' @examples
#' is_kernel(gaussian)
#' no_kernel <- function(x) 1
#' is_kernel(no_kernel)
#'
#' @include integrable_function.R
#' @export
validate_Kernel <- function(obj){
  stopifnot("object has to be of class Kernel"=inherits(obj, "Kernel"))
  validate_IntegrableFunction(obj)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]
  stopifnot("The integral of a kernel over its support has to be one"=
              (abs(integrate(object, lower = lower, upper = upper)[[1]] - 1)
         < integrate(object, lower = lower, upper = upper)[[2]]))
  invisible(obj)
}

#' @include integrable_function.R
new_Kernel <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}
