#' Shifting kernels to a coefficient and applying a bandwidth.
#'
#' @description
#' The \code{kernel_transform} is shifting and applying a bandwidth on a kernel function.
#'
#' @param kernel a kernel object.
#' @param sample the observation/coefficient as a numerical scalar.
#' @param h the bandwidth as a numerical scalar.
#' @return the transformed, i.e. shifted and scaled kernel function.
#'
#' @examples
#' # shifting the gaussian kernel
#' shifted_gaussian <- kernelTransform(gaussian, 2, 1)
#' # stretching the gaussian kernel
#' stretched_gaussian <- kernelTransform(gaussian, 0, 2)
#' x <- seq(from = -5, to = 5, length.out = 1000)
#'
#' plot(x, gaussian$fun(x),
#'      xlim=c(-5,5), ylim=c(0,1),
#'      main="Kernels", xlab="", ylab="",
#'      col="black", type="l")
#' lines(x, shifted_gaussian(x), col="red")
#' lines(x, stretched_gaussian(x), col="blue")
#' legend("topright",
#'        legend=c("gaussian", "shifted_gaussian","stretched_gaussian"),
#'        col=c("black","red", "blue"), lty=1, cex=0.8)
#'
#' @export
kernel_transform <- function(kernel, sample, h) {
  # Kernel conditions
  tryCatch({validate_Kernel(kernel)}, error="the kernel has to be valid")

  # x_i conditions
  stopifnot(is.numeric(sample))
  stopifnot(length(sample) == 1)

  # Bandwidth conditions
  stopifnot(is.numeric(h))
  stopifnot(length(h) == 1)
  stopifnot(h > 0)

  force(kernel)
  force(sample)
  force(h)

  fun <- function(x){
    kernel$fun((x-sample)/h)/h
  }
  support <- (h * kernel$support + sample)

  Kernel(fun, support)
}

