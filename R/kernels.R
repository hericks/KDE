#' Kernel functions.
#'
#' @description
#' Built-in kernel functions.
#'
#' @details{
#' A kernel is a real valued, integrable function, such that its integral over the real numbers equals one.
#' The built-in kernel functions are vectorised.
#' \describe{
#' \strong{List of built-in kernels:}
#'   \item{\code{rectangular}}
#'   \item{\code{triangular}}
#'   \item{\code{epanechnikov}}
#'   \item{\code{biweight}}
#'   \item{\code{triweight}}
#'   \item{\code{tricube}}
#'   \item{\code{gaussian}}
#'   \item{\code{cosine}}
#'   \item{\code{logistic}}
#'   \item{\code{sigmoidFunction}}
#'   \item{\code{silverman}}
#'   }
#' }
#' @param u vector containing numerical values.
#'
#' @return returning evaluation of the kernel function in u.
#'
#' @examples
#' x <- seq(from = -4, to = 4, length.out = 1000)
#' plot(x, gaussian(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x,triangular(x), col="red")
#' lines(x, rectangular(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "triangular","rectangular"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#' @name kernels
#' @export
NULL


checkKernelConditions <- function(u){
  stopifnot("the argument of a kernel has to be numeric"=is.numeric(u))
}

rectangular <- function(u){
  checkKernelConditions(u)
  return(1/2*(abs(u) <= 1))
}

triangular <- function(u){
  checkKernelConditions(u)
  return((1 - abs(u)) * (abs(u) <= 1))
}

epanechnikov <- function(u){
  checkKernelConditions(u)
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

biweight <- function(u){
  checkKernelConditions(u)
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

triweight <- function(u){
  checkKernelConditions(u)
  return((35/32 * (1 - u^2)^3) * (abs(u) <= 1))
}

tricube <- function(u){
  checkKernelConditions(u)
  return(70/81 * (1 - abs(u^3))^3 * (abs(u) <= 1))
}

gaussian <- function(u){
  checkKernelConditions(u)
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

cosine <- function(u){
  checkKernelConditions(u)
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

logistic <- function(u){
  checkKernelConditions(u)
  return(1/(exp(u) + 2 + exp(-u)) )
}

sigmoidFunction <- function(u){
  checkKernelConditions(u)
  return(2/pi * (1/(exp(u) + exp(-u)) ) )
}

silverman <- function(u){
  checkKernelConditions(u)
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}



