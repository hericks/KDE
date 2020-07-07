#' Validation if given function is probability density function.
#'
#' @description The \code{is_density()} function is used to check if a given
#'   function satisfies the conditions of a probability density function. It has
#'   to be a integrable, nonnegativ numerical function with integral over the
#'   real numbers equal to 1.
#'
#' @param d the fuction to be tested.
#'
#' @return Boolean. Returns True if d is a probability density function, FALSE if not.
#'   Is not supossed to return an error.
#'
#' @examples
#' TODO
#'

is_density <- function(d) {
  if (class(d) != "function") return(FALSE)

  #d is nonnegative
  neg_d <- Vectorize(function(x) max(0, -d(x)))
  if(isFALSE(tryCatch(integrate(neg_d, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  neg_integral_d <- integrate(neg_d, -Inf, Inf)[[1]]
  if (neg_integral_d != 0) return(FALSE)

  # d is normalized
  pos_d <- Vectorize(function(x) max(0, d(x)))
  if(isFALSE(tryCatch(integrate(pos_d, -Inf, Inf), error=function(e) FALSE))){
    return(FALSE)
  }
  pos_integral_d <- integrate(pos_d, -Inf, Inf)[[1]]

  if (pos_integral_d != 1) return(FALSE)
}

