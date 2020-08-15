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
