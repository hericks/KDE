#' Shifting kernels to a coefficient and applying a bandwidth.
#'
#' @description
#' The \code{kernel_transform} is shifting and applying a bandwidth to a kernel.
#'
#' @param kernel a kernel object.
#' @param sample the observation/coefficient as a numerical scalar.
#' @param h the bandwidth as a numerical scalar.
#' @param subdivisions a integer vector of length 1, used for the subdivisions parameter for the function \code{\link[stats:integrate]{integrate}} in \code{\link[KDE:Kernel]{Kernel}}.
#'
#' @return The transformed (shifted and scaled) kernel as valid S3 object of
#'   class \code{Kernel}.
#'
#' @seealso
#' \code{\link[KDE:Kernel]{Kernel}}
#' \code{\link[KDE:validate_Kernel]{validate_Kernel}}
#' \code{\link[KDE:IntegrableFunction]{IntegrableFunction}}
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
kernel_transform <- function(kernel, sample, h, subdivisions=100L) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # x_i conditions
  stopifnot(is.numeric(sample))
  stopifnot(length(sample) == 1)

  # Bandwidth conditions
  stopifnot(is.numeric(h))
  stopifnot(length(h) == 1)
  stopifnot(h > 0)

  Kernel(function(x) kernel$fun((x-sample)/h)/h, h*kernel$support + sample, subdivisions=subdivisions)
}

