#' Validation if object is a kernel function.
#'
#' @description
#' The \code{is_kernel()} function is used for validating a \link[KDE:kernels]{kernel function}.
#' It has to be a integrable, numerical function with integral over R equal to one.
#'
#' @param object the fuction to be tested.
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
#' @export
is_kernel <- function(object) {
  if (class(object) != "function") return(FALSE)

  # neg/pos part of the kernel function to be tested
  pos_object <- Vectorize(function(x) max(0, object(x)))
  neg_object <- Vectorize(function(x) max(0, -object(x)))

  # is_kernel returns a boolean and is not supposed to throw an error
  if(isFALSE(tryCatch(integrate(pos_object, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  if(isFALSE(tryCatch(integrate(neg_object, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  pos_integral <- integrate(pos_object, -Inf, Inf)[[1]]
  neg_integral <- integrate(neg_object, -Inf, Inf)[[1]]

  if(is.infinite(pos_integral) && is.infinite(pos_integral)) return(FALSE)

  isTRUE(abs(integrate(object, -Inf, Inf)[[1]] - 1) < integrate(object, -Inf, Inf)[[2]])
}




