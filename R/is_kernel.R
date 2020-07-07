#' Validation if object is a kernel function.
#'
#' @description
#' The \code{is_kernel()} function is used for validating a \link[KDE:kernels]{kernel function}.
#' It has to be a integrable, numerical function with integral over R equal to one.
#'
#' @details
#' The is_kernel function works with centered kernels, searching a non-zero interval around the center,
#' checking integrability over that interval.
#' The function is checking the points
#' (-1e1, -1e0, -1e-1, -1e-2, -1e-3, -1e-4, -1e-5, -1e-6, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1) around the center.
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
#' @export
is_kernel <- function(object, center=0) {
  #TODO Argchecks
  if(class(object) != "function") return(FALSE)
  if(!is.numeric(object(center))) return(FALSE)

  force(object)
  force(center)
  if(isFALSE(tryCatch({find_borders(object,center)}, error=function(e) FALSE))){
    return(FALSE)
  }
  borders <- find_borders(object, center)
  if(!(is_integrable(object, borders[1], borders[2]))) return(FALSE)

  isTRUE(abs(integrate(object, borders[1], borders[2])[[1]] - 1) < integrate(object, borders[1], borders[2])[[2]])
}

is_integrable <- function(object, lower=-Inf, upper=Inf){
  # neg/pos part of the kernel function to be tested
  pos_object <- Vectorize(function(x) max(0, object(x)))
  neg_object <- Vectorize(function(x) max(0, -object(x)))

  # is_kernel returns a boolean and is not supposed to throw an error
  if(isFALSE(tryCatch({
    pos_int <- integrate(pos_object, lower=lower, upper=upper)
    neg_int <- integrate(neg_object, lower=lower, upper=upper)
  }, error=function(e) FALSE))) {
    return(FALSE)
  }
  pos_integral <- integrate(pos_object, lower, upper)[[1]]
  neg_integral <- integrate(neg_object, lower, upper)[[1]]
  if(is.infinite(pos_integral) && is.infinite(pos_integral)){ return(FALSE)}

  return(TRUE)
}

find_borders <- function(fun, center=0){
  x <- c(1e1, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
  y <- -x
  x <- x + center
  y <- y + center

  check_x <- fun(x) != 0
  check_y <- fun(y) != 0
  if(!(any(check_x) && any(check_y))){
    return(c(-1e-6,1e-6))
  }
  non_zero <- min(c(which(check_x==TRUE), which(check_y==TRUE)))
  if(non_zero != 1){
    return(c(y[non_zero-1], x[non_zero-1]))
  }
  else{
    return(c(-Inf,Inf))
  }
}
