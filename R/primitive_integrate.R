integrate_primitive <- function(integrand,
                                lower,
                                upper,
                                stepsize,
                                max_length) {
  # TODO: Conditions

  max_length <- ceiling(max_length)
  integration_length <- upper - lower
  num_steps <- ceiling(integration_length/stepsize + 1)

  if (num_steps <= max_length) {
    eval_points <- seq(from=lower, to=upper, by=stepsize)
  } else {
    eval_points <- seq(from=lower, to=upper, length.out = max_length)
    used_step_size <- ifelse(max_length == 1, Inf, eval_points[2] - eval_points[1])
    warning(paste0("Warning (primitive_integrate): Had to increase step size from ", stepsize, " to ", used_step_size, "."), call. = FALSE)
  }

  eval_values <- integrand(eval_points)
  sum(eval_values)*integration_length/length(eval_points)
}
