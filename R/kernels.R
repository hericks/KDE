#' Kernel functions
#'
#' @description
#' Built-in kernel functions.
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
NULL

checkKernelConditions <- function(u){
  stopifnot("the argument of a kernel has to be numeric"=is.numeric(u))
}

#' Rectangular Function
#' @description
#' The rectangular function.
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
rectangular <- function(u){
  checkKernelConditions(u)
  return(1/2*(abs(u) <= 1))
}

#' Triangular Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
triangular <- function(u){
  checkKernelConditions(u)
  return((1 - abs(u)) * (abs(u) <= 1))
}

#' Epanechnikov Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
epanechnikov <- function(u){
  checkKernelConditions(u)
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

#' Biweight Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
biweight <- function(u){
  checkKernelConditions(u)
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

#' Triweight Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
triweight <- function(u){
  checkKernelConditions(u)
  return((35/32 * (1 - u^2)^3) * (abs(u) <= 1))
}

#' Tricube Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
tricube <- function(u){
  checkKernelConditions(u)
  return(70/81 * (1 - abs(u^3))^3 * (abs(u) <= 1))
}

#' Gaussian Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
gaussian <- function(u){
  checkKernelConditions(u)
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

#' Cosine Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
cosine <- function(u){
  checkKernelConditions(u)
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

#' Logistic Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
logistic <- function(u){
  checkKernelConditions(u)
  return(1/(exp(u) + 2 + exp(-u)) )
}

#' Sigmoid Function
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
sigmoidFunction <- function(u){
  checkKernelConditions(u)
  return(2/pi * (1/(exp(u) + exp(-u)) ) )
}

#' Silverman Function
#'
#' @param u vector of numerical values
#'
#' @family kernels
#' @seealso
#' \code{\link[KDE:kernels]{kernels}} for more information about kernels.
#' @export
silverman <- function(u){
  checkKernelConditions(u)
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}



