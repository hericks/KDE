check_kernel_conditions <- function(u){
  stopifnot("the argument of a kernel has to be numeric"=is.numeric(u))
}

rectangular_function <- function(u){
  check_kernel_conditions(u)
  return(1/2*(abs(u) <= 1))
}

#' Rectangular Kernel
#'
#' @description The rectangular kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the rectangular function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
rectangular <- Kernel(rectangular_function, c(-1,1))

triangular_function <- function(u){
  check_kernel_conditions(u)
  return((1 - abs(u)) * (abs(u) <= 1))
}

#' Triangular Kernel
#'
#' @description The triangular kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the triangular function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link[KDE:Kernel]{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
triangular <- Kernel(triangular_function, c(-1,1))


epanechnikov_function <- function(u){
  check_kernel_conditions(u)
  return(3/4 * (1 - u^2) * (abs(u) <= 1))
}

#' Epanechnikov Kernel
#'
#' @description The epanechnikov kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the epanechnikov function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
epanechnikov <- Kernel(epanechnikov_function, c(-1,1))

biweight_function <- function(u){
  check_kernel_conditions(u)
  return(15/16 * (1 - u^2)^2 * (abs(u) <= 1))
}

#' Biweight Kernel
#'
#' @description The biweight kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the biweight function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso
#' \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
biweight <- Kernel(biweight_function, c(-1,1))

triweight_function <- function(u){
  check_kernel_conditions(u)
  return((35/32 * (1 - u^2)^3) * (abs(u) <= 1))
}

#' Triweight Kernel
#'
#' @description The triweight kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the triweight function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
triweight <- Kernel(triweight_function, c(-1,1))

tricube_function <- function(u){
  check_kernel_conditions(u)
  return(70/81 * (1 - abs(u^3))^3 * (abs(u) <= 1))
}

#' Tricube Kernel
#'
#' @description The tricube kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the tricube function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
tricube <- Kernel(tricube_function, c(-1,1))

gaussian_function <- function(u){
  check_kernel_conditions(u)
  return(1/sqrt(2 * pi) * exp(-1/2 * u^2))
}

#' Gaussian Kernel
#'
#' @description The gaussian kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the gaussian function
#'
#' * \bold{\code{support}}: \code{c(-10,10)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
gaussian <- Kernel(gaussian_function, c(-10,10))

cosine_function <- function(u){
  check_kernel_conditions(u)
  return(pi/4 * cos(pi/2 * u) * (abs(u) <= 1))
}

#' Cosine Kernel
#'
#' @description The cosine kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the cosine function
#'
#' * \bold{\code{support}}: \code{c(-1,1)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
cosine <- Kernel(cosine_function, c(-1,1))

logistic_function <- function(u){
  check_kernel_conditions(u)
  return(1/(exp(u) + 2 + exp(-u)) )
}

#' Logistic Kernel
#'
#' @description The logistic kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the logistic function
#'
#' * \bold{\code{support}}: \code{c(-20,20)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
logistic <- Kernel(logistic_function, c(-20, 20))

sigmoid_function <- function(u){
  check_kernel_conditions(u)
  return(2/pi * (1/(exp(u) + exp(-u))))
}

#' Sigmoid Kernel
#'
#' @description The sigmoid kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the sigmoid function
#'
#' * \bold{\code{support}}: \code{c(-20,20)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
sigmoid <- Kernel(sigmoid_function, c(-20,20))

silverman_function <- function(u){
  check_kernel_conditions(u)
  return(0.5 * exp(-abs(u)/sqrt(2)) * sin(abs(u)/sqrt(2) + pi/4))
}

#' Silverman Kernel
#'
#' @description The silverman kernel is S3 object of class \code{Kernel}
#'   provided by the \code{KDE} package.
#'
#' @format
#' An object of S3 class \code{\link{Kernel}} with named entries
#'
#' * \bold{\code{fun}}: the silverman function
#'
#' * \bold{\code{support}}: \code{c(-25,25)}
#'
#' * \bold{\code{subdivisions}}: \code{1000L}.
#'
#' @family kernels
#'
#' @seealso \code{\link{Kernel}} for more information about kernels.
#'
#' @include kernel.R
#' @export
silverman <- Kernel(silverman_function, c(-25,25))
