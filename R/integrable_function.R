#' Integrable-Function Classes
#'
#' @description
#' The S3 class \code{IntegrableFunction} serves as superclass for the S3 classes \code{Density} and \code{Kernel}.
#'
#' @export
IntegrableFunction <- function(fun, support){
  func <- new_IntegrableFunction(fun, support)
  validate_IntegrableFunction(func)
  func
}


#' TODO: documentation
#'
#'@export
validate_IntegrableFunction <- function(x){
  # to test the basic structure
  stopifnot("Object must inherit 'IntegrableFunction'"=inherits(x, "IntegrableFunction"))
  stopifnot("Object must be a list"=is.list(x))
  stopifnot("Object must contain entry 'fun'"="fun" %in% names(x))
  stopifnot("Object must contain entry 'support'"="support" %in% names(x))

  # to test the structure of 'fun' entry
  fun <- x$fun
  stopifnot("Entry 'fun' must be a function"=is.function(fun))

  # to test the structure of 'support' entry
  support <- x$support
  stopifnot("Entry 'support' must be numeric"=is.numeric(support))
  stopifnot("Entry 'support' must be of length 2"=identical(length(support), 2L))
  stopifnot("Entry 'support' must contain lower- before upperbound"=support[1] < support[2])

  # TODO: test that fun is equal to zero outside of support

  # to test the the integrability on the given support
  tryCatch({
    integration_value <- integrate(function(x) abs(fun(x)), support[1], support[2])[[1]]
  }, error=function(cond) {
    stop(paste("Failed to integrate the function:", cond))
  })

  stopifnot("The integral of the absolute function must integrate to a finite value"=is.finite(integration_value))
  invisible(x)
}

# This function is only called by the public IntegrableFunction constructor followed by a call of validate_IntegrableFunction.
# Further checking is done by validate_IntegrableFunction.
new_IntegrableFunction <- function(fun, support=NULL, ..., subclass=NULL){
  stopifnot("class of fun has to be function"=is.function(fun))

  if(is.null(support)){
    support <- find_support(fun)
  }

  structure(
    list("fun"=fun, "support"=support),
    ...,
    class=c(subclass, "IntegrableFunction")
  )
}

find_support <- function(fun) {
  testing_points <- c(-10**(10:-10), 10**(-10:10))
  testing_values <- fun(testing_points)

  stopifnot("fun has to return numerical values"=is.numeric(testing_values))
  stopifnot("fun hat to be vectorised"=identical(length(testing_values), length(testing_points)))

  # can't use isFALSE, since all.equal return value is TRUE or a vector of mode "character"
  non_zero_indices <- which(sapply(testing_values, function(x) !isTRUE(all.equal(x, 0))))

  if (length(non_zero_indices) == 0L) return(c(-Inf, Inf))

  lower_index <- min(non_zero_indices) - 1
  upper_index <- max(non_zero_indices) + 1

  lower_bound <- ifelse(lower_index < 1, -Inf, testing_points[lower_index])
  upper_bound <- ifelse(upper_index > length(testing_points), Inf, testing_points[upper_index])

  c(lower_bound, upper_bound)
}
