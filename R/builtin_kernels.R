check_kernel_conditions <- function(u){
  stopifnot("the argument of a kernel has to be numeric"=is.numeric(u))
}

rectangular_function <- function(u){
  check_kernel_conditions(u)
  return(1/2*(abs(u) <= 1))
}

#' Rectangular Function
#'
#' @description
#' The rectangular function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
rectangular <- Kernel(rectangular_function, c(-Inf,Inf))

triangular_function <- function(u){
  check_kernel_conditions(u)
  return((1 - abs(u)) * (abs(u) <= 1))
}

#' Triangular Function
#' @description
#' The triangular function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
triangular <- Kernel(triangular_function, c(-Inf,Inf))


epanechnikov_function <- function(u){
  check_kernel_conditions(u)
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

#' Epanechnikov Function
#'
#' @description
#' The epanechnikov function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
epanechnikov <- Kernel(epanechnikov_function, c(-Inf,Inf))

biweight_function <- function(u){
  check_kernel_conditions(u)
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

#' Biweight Function
#'
#' @description
#' The biweight function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
biweight <- Kernel(biweight_function, c(-Inf,Inf))

triweight_function <- function(u){
  check_kernel_conditions(u)
  return((35/32 * (1 - u^2)^3) * (abs(u) <= 1))
}

#' Triweight Function
#'
#' @description
#' The triweight function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
triweight <- Kernel(triweight_function, c(-Inf,Inf))

tricube_function <- function(u){
  check_kernel_conditions(u)
  return(70/81 * (1 - abs(u^3))^3 * (abs(u) <= 1))
}

#' Tricube Function
#'
#' @description
#' The tricube function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
tricube <- Kernel(tricube_function, c(-Inf,Inf))

gaussian_function <- function(u){
  check_kernel_conditions(u)
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

#' Gaussian Function
#'
#' @description
#' The gaussian function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
gaussian <- Kernel(gaussian_function, c(-Inf,Inf))

cosine_function <- function(u){
  check_kernel_conditions(u)
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

#' Cosine Function
#'
#' @description
#' The cosine function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
cosine <- Kernel(cosine_function, c(-Inf,Inf))

logistic_function <- function(u){
  check_kernel_conditions(u)
  return(1/(exp(u) + 2 + exp(-u)) )
}

#' Logistic Function
#'
#' @description
#' The logistic function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
logistic <- Kernel(logistic_function, c(-Inf,Inf))

sigmoid_function <- function(u){
  check_kernel_conditions(u)
  return(2/pi * (1/(exp(u) + exp(-u)) ) )
}

#' Sigmoid Function
#'
#' @description
#' The sigmoid function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
sigmoidFunction <- Kernel(sigmoid_function, c(-Inf,Inf))

silverman_function <- function(u){
  check_kernel_conditions(u)
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}

#' Silverman Function
#'
#' @description
#' The silverman function is used as a \link[KDE:Kernel]{kernel}.
#'
#' @param u vector of numerical values
#'
#' @family kernels
#'
#' @seealso
#' \code{\link[KDE:Kernel]{kernels}} for more information about kernels.
#'
#' @include kernel.R
#' @export
silverman <- Kernel(silverman_function, c(-Inf,Inf))





