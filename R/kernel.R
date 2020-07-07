#' Kernel
#'
#' @description
#' Built-in kernels are S3 Objects, holding a function and its support.
#'
#' @details{
#' A kernel is a real valued, integrable function, such that its integral over the real numbers equals one.
#' The built-in kernel functions are vectorised.
#' \describe{
#' \strong{List of built-in kernels:}
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
#' @param u vector containing numerical values.
#'
#' @return returning evaluation of the kernel function in u.
#'
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
#'
#' @export
Kernel <- function(fun, support) {
  kern <- new_Kernel(fun, support)
  validate_Kernel(kern)
  kern
}

#' Validation if object is a kernel function.
#'
#' @description
#' The \code{is_kernel()} function is used for validating a \link[KDE:Kernel]{Kernel}.
#' It has to be a integrable, numerical function with integral over R equal to one.
#'
#' @details
#' The is_kernel function works with centered kernels, searching a non-zero interval around the center,
#' checking integrability over that interval.
#'
#' @param object the object to be tested.
#' @param center the center of the function object. Default value is zero.
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
#' @include evaluate.R
#' @export
validate_Kernel <- function(obj){
  if(!inherits(obj, "Kernel")) return(FALSE)
  validate_IntegrableFunction(obj)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]
  isTRUE(abs(integrate(evaluate_safe(obj), lower = lower, upper = upper)[[1]] - 1)
         < integrate(evaluate_safe(obj), lower = lower, upper = upper)[[2]])
}

new_Kernel <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Kernel"))
}
