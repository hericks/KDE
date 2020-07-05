#' Validation if kernel is a kernel function.
#'
#' @description
#' The \code{is_kernel()} function is used for validating a \link[KDE:kernels]{kernel function}.
#' It has to be a integrable, numerical function with integral over R equal to one.
#'
#' @param kernel the fuction to be tested.
#'
#' @return Boolean. Returns True if kernel is a kernel function, FALSE if not. Is not supossed to return an error.
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
#' @export
is_kernel <- function(kernel) {
  if (class(kernel) != "function") return(FALSE)

  pos_kernel <- Vectorize(function(x) max(0, kernel(x)))
  neg_kernel <- Vectorize(function(x) max(0, -kernel(x)))

  # is_kernel returns a boolean and is not supposed to throw an error
  if(isFALSE(tryCatch(integrate(pos_kernel, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  if(isFALSE(tryCatch(integrate(neg_kernel, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  pos_integral <- integrate(pos_kernel, -Inf, Inf)[[1]]
  neg_integral <- integrate(neg_kernel, -Inf, Inf)[[1]]

  if(is.infinite(pos_integral) && is.infinite(pos_integral)) return(FALSE)

  isTRUE(abs(integrate(kernel, -Inf, Inf)[[1]] - 1) < integrate(kernel, -Inf, Inf)[[2]])
}




