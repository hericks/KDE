#' Shifting Kernels to a Coefficient and applying a Bandwidth.
#'
#' @description
#' \code{kernel_transform} is shifting and applying a bandwidth to a kernel.
#'
#' @param kernel a kernel as S3 object of the class \link{Kernel}.
#' @param sample numeric scalar; the observation.
#' @param bandwidth numeric scalar; the bandwidth.
#' @param subdivisions positive numeric scalar; subdivisions parameter
#'   internally passed to \code{\link{integrate_primitive}}.
#'
#' @details The validation of the returned estimator as
#'   \code{\link{IntegrableFunction}} relies on the function
#'   \code{\link{integrate_primitive}}, thus the \code{subdivisions} parameter.
#'
#' @return The transformed (shifted and scaled) kernel as valid S3 object of
#'   class \code{Kernel}.
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @examples
#' # shifting the gaussian kernel
#' shifted_gaussian <- kernel_transform(gaussian, 2, 1)
#' # stretching the gaussian kernel
#' stretched_gaussian <- kernel_transform(gaussian, 0, 2)
#' x <- seq(from = -5, to = 5, length.out = 1000)
#'
#' plot(x, gaussian$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x, shifted_gaussian$fun(x), col="red")
#' lines(x, stretched_gaussian$fun(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "shifted_gaussian","stretched_gaussian"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#' @export
kernel_transform <- function(kernel, sample, bandwidth, subdivisions=1000L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # x_i conditions
  stopifnot(is.numeric(sample))
  stopifnot(length(sample) == 1)

  # Bandwidth conditions
  stopifnot(is.numeric(bandwidth))
  stopifnot(length(bandwidth) == 1)
  stopifnot(bandwidth > 0)

  Kernel(function(x) kernel$fun((x-sample)/bandwidth)/bandwidth, bandwidth*kernel$support + sample, subdivisions=subdivisions)
}

