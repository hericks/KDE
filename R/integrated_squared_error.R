#' @include integrate_primitive.R
#'
#' @export
integrated_squared_error <- function(f1, f2, support, subdivisions= 1000L) {
  integrate_primitive(integrand = function(x) {(f1(x) - f2(x))^2},
                      lower = support[1],
                      upper = support[2],
                      subdivisions = subdivisions
  )
}

#' @include integrate_primitive.R
#'
#' @export
integrated_squared_error_kde <- function(density, kernel, samples, bandwidth, subdivisions) {
  estimator <- kernel_density_estimator(kernel, samples, bandwidth, subdivisions)
  integrated_squared_error(estimator$fun,
                           density$fun,
                           support=range(estimator$support, density$support),
                           subdivisions = subdivisions)$value
}
