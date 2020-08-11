#' Construct a kernel density estimator
#'
#' @description The \code{kernel_density_estimator} function that approximates a probability density function from given sampels.
#'
#' @param kernel A S3 object of the class \link[KDE:Kernel]{Kernel}.
#' @param samples A numerical vector of observations.
#' @param bandwidth A non-negative numeric value to use as the bandwidth for the
#'   kernel.
#' @param subdivisions A integer vector of length 1 used for the subdivisions parameter of the builtin R-function \code{\link[stats:integrate]{integrate}}.
#'
#' @details{
#' * \strong{`kernel`} Kernels in this package are S3 objects of the class \link[KDE:Kernel]{Kernel}.\cr
#' See \link[KDE:Kernel]{Kernel} for more informations on kernels.
#'
#' * \strong{`subdivisions`} is a integer value used for the subdivisions parameter for \code{\link[stats:integrate]{integrate}}.
#' The subdivisions parameter is required to be large enough, such that \code{\link[stats:integrate]{integrate}} can work properly.
#' The default value is set to 100L. Be aware that too large numbers can cause long runtimes!
#'
#' For more information about kernel density estimaton, see "Nonparametric Estimation" by Fabienne Comte.
#' }
#'
#' @return A S3 object of class \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}.
#'
#' @source
#' \itemize{Nonparametric Estimation, Comte `[`2017`]`, ISBN: 978-2-36693-030-6}
#'
#' @seealso
#' \code{\link[KDE:Kernel]{Kernel}}
#' \code{\link[KDE:validate_Kernel]{validate_Kernel}}
#' \code{\link[KDE:Kernel]{IntegrableFunction}}
#'
#' @export
kernel_density_estimator <- function(kernel, samples, bandwidth=1, subdivisions=100L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # Samples conditions
  stopifnot(is.numeric(samples))
  stopifnot(length(samples) > 0)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  # Subdivisions conditions
  subdivisions <- subdivisions
  stopifnot("Entry 'subdivisions' must be numeric"=is.numeric(subdivisions))
  stopifnot("Entry 'subdivisions' must be positive"=subdivisions > 0)

  kernel_eval <- kernel$fun

  estimator_eval <- function(x) {
    stopifnot("input has to be numeric"=is.numeric(x))

    ret <- numeric(length(x))
    for (x0 in samples) {
      ret <- ret + kernel_eval((x - x0)/bandwidth)/bandwidth
    }
    ret/length(samples)
  }

  support <- c(bandwidth*kernel$support[1] + min(samples), bandwidth*kernel$support[2] + max(samples))

  IntegrableFunction(estimator_eval, support, subdivisions=subdivisions)
}







