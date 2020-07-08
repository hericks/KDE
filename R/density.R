#'
#' @include integrable_function.R
#'
#' @export
Density <- function(fun, support = NULL){
  den <- new_Density(fun, support)
  validate_Density(den)
  den
}

#'
#'@include integrable_function.R
#'
#' @export
validate_Density <- function(obj){
  stopifnot("object has to be of class Density"=inherits(obj, "Density"))
  validate_IntegrableFunction(obj)
  object <- obj$fun
  lower <- obj$support[1]
  upper <- obj$support[2]

  testing_points <- seq(max(lower, -1e10), min(upper, 1e10), length.out=1e4)
  # prevent double checking with validate_integrable function
  stopifnot("density functions are non-negative" = all(object(testing_points) >= 0))

  stopifnot("The integral of a density over its support has to be one"=
              isTRUE(all.equal(integrate(object, lower = lower, upper = upper)[[1]], 1)))
  invisible(obj)
}

new_Density <- function(fun, support, ..., subclass=NULL){

  new_IntegrableFunction(fun,support, subclass=c(subclass,"Density"))
}
