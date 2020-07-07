#' Integrable-Function Classes
#'
#' @description The S3 class \code{IntegrableFunction} serves as superclass for
#' the S3 classes \code{Density} and \code{Kernel}. The constructor tries to
#' construct a valid \code{IntegrableFunction}.
#'
#' @param fun an \code{R} function taking a single numeric argument and
#'   returning a numeric vector of the same length: see 'Details'.
#' @param support a numerical vector of length 2 containing the lower- and
#'   upperbound in the first and second entry respectively.
#'   \code{IntegrableFunction} will try to find bounds on the support itself if
#'   \code{NULL} is passed: see 'Details'.
#'
#' @details In this package the S3 class \code{IntegrableFunction} exists to ensure
#' the integrability of \code{R} functions. Valid \code{IntegrableFunction} objects are build
#' on lists and contain two names entries:
#'
#' * \strong{`fun`} is an \code{R} function taking a single numeric argument and
#' returning a numeric vector of the same length. This function should return
#' zero outside of the interval given in the \code{support} entry.
#'
#' * \strong{`support`} is a numeric vector of length 2 containing a lower- and
#' upperbound for the support of the function stored in \code{fun} in its first
#' and second entry respectively. In particular the values \code{-Inf} and
#' \code{Inf} are allowed.
#'
#' The constructor \code{IntegrableFunction} tries to construct an
#' \code{IntegrableFunction} object based on the provided arguments. In
#' particular it assigns the correct class attribute and verifies the correct
#' structure of the \code{fun} and \code{support} entries, such that the returned
#' object passes a test using [validate_IntegrableFunction()].
#'
#' @export
IntegrableFunction <- function(fun, support=NULL){
  func <- new_IntegrableFunction(fun, support)
  validate_IntegrableFunction(func)
  func
}


#' Validate Integrable-Function Class
#'
#' @description This function serves as a validator for the S3 class \code{IntegrableFunction}.
#' See 'Details' for further information and potential flaws.
#'
#' @param x The object to validate as S3 object of class \code{IntegrableFunction}.
#'
#' @details This function tries to verify an object as S3 object of class
#'   [IntegrableFunction]. In particular the formal structure is ensured.
#'   Additionally this function tries to validate the additional conditions of
#'   valid \code{IntegrableFunction} objects, such that
#'
#'   1. The \code{fun} entry is vectorised in its argument, returning numerical
#'   values only.
#'
#'   2. The \code{fun} entry returns zeros for inputs outside of
#'   the support interval specified in the \code{support} entry.
#'
#'   3. \code{integrate} can be used to integrate \code{fun} over the interval
#'   \code{entry} without throwing an error.
#'
#'  Like all numerical routines, this function can evaluate represented function
#'  on a finite set of points only. If the represented function returns valid
#'  results over nearly all its range, it is possible that this function misses
#'  unexpected/wrong return values. Thus, using [IntegrableFunction] or
#'  [validate_IntegrableFunction] to construct and validate objects representing
#'  integrable functions is \emph{not} sufficient to ensure the listed
#'  properties \[1-3\], but more of a sanity-check.
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
