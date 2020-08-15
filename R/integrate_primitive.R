#' Primitive Integration
#'
#' @description \code{integrate_primitive} integrates an in the mathematical
#'   sense integrable function over a compact interval using a basic
#'   step-function approximation.
#'
#' @param integrand a \code{R} function taking a single numeric argument and
#'   returning a numeric vector of the same length. Further the function should
#'   be integrable in the mathematical sense.
#' @param lower the lower bound for the integration. Must be finite.
#' @param upper the upper bound for the integration. Must be finite.
#' @param subdivisions the number of evaluation points used in the procedure.
#' @param check a boolean value to activate the collection of an additional measure of convergence.
#'
#' @details The \code{KDE} package has to frequently handle functions, which are
#'   integrable in the mathematical sense but are hard to integrate consistently
#'   using the base \code{R} function \code{integrate} and equivalent choices.
#'   The typical integration methods often fail for non-smooth functions or functions
#'   with a support consisting of multiple narrow disjoint intervals.
#'
#'   Since the function to integrate is often computationally cheap to evaluate,
#'   approximating the function on a narrow grid and integrating the approximation
#'   as primitive function yields a practical approach.
#'
#'   \code{integrate_primitive} makes use of this approach. The support has to
#'   be finite to allow for a discretisation grid. In particular \code{lower}
#'   and \code{upper} must be finite. The number of evaluation points is
#'   specified in the positive numerical parameter \code{subdivisions}. In
#'   particular a higher number of subdivisions yields a more accurate result at
#'   the expense of longer runtimes.
#'
#'   The \code{check} parameter can be set to enable the collection of an
#'   additional measure of convergence \code{relError}. The relative error is
#'   calculated by computing the integral for increasing numbers of subdivisions
#'   and comparing the successive results.
#'
#' @return \code{integrate_primitive} returns a named list containing the entries
#'   \code{value} and \code{relError} for the integration value and relative
#'   Error respectively. If the \code{check} parameter is not set
#'   \code{relError} will always be equal to \code{NULL}.
#'
#' @export
integrate_primitive <- function(integrand,
                                lower,
                                upper,
                                subdivisions = 1000L,
                                check = FALSE) {
  # TODO: Conditions
  integration_length <- upper - lower

  eval_points <- seq(from=lower, to=upper, length.out = subdivisions)
  eval_values <- integrand(eval_points)
  eval_values <- eval_values[!is.infinite(eval_values)]

  integral <- sum(eval_values)*integration_length/length(eval_points)
  rel_error <- NULL

  if (check) {
    eval_points <- seq(from=lower, to=upper, length.out = ceiling(subdivisions/2))
    eval_values <- integrand(eval_points)
    eval_values <- eval_values[!is.infinite(eval_values)]

    approx_integral <- sum(eval_values)*integration_length/length(eval_points)

    abs_error <- abs(integral - approx_integral)
    rel_error <- ifelse(all.equal(abs_error, 0), 0, abs_error/abs(integral))
  }

  return(list(value=integral, relError=rel_error))
}