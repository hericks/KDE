#' Kernel functions.
#'
#' @description
#' Listing of the built-in kernel functions. A kernel is a real valued function \eqn{K: \mathbb{R} \rightarrow \mathbb{R}}, integrable and such that \eqn{\int\limits_{-\infty}^{\infty} K(u) du = 1}.
#' The built-in kernel functions are vectorised.
#'
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
#' @export
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



