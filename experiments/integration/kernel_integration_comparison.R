kernels <- list(rectangular,
                triangular,
                epanechnikov,
                biweight,
                triweight,
                gaussian,
                cosine,
                logistic,
                sigmoid,
                silverman)

# Settings
subdivisions = 1000L

# Error Analysis
primitive_abs_errors <- double(length(kernels))
for (i in seq_along(kernels)) {
  kernel <- kernels[[i]]
  I <- integrate_primitive(integrand = kernel$fun,
                           lower = kernel$support[1],
                           upper = kernel$support[2],
                           subdivisions = subdivisions)$value

  primitive_abs_errors[i] <- abs(I - 1)
}

base_abs_errors <- double(length(kernels))
for (i in seq_along(kernels)) {
  kernel <- kernels[[i]]
  I <- integrate(f = kernel$fun,
                 lower = kernel$support[1],
                 upper = kernel$support[2])$value

  base_abs_errors[i] <- abs(I - 1)
}

summary(primitive_abs_errors)
summary(base_abs_errors)

# Performance Analysis
primitive_integration <- function() {
  for (kernel in kernels) {
    integrate_primitive(integrand = kernel$fun,
                        lower = kernel$support[1],
                        upper = kernel$support[2],
                        subdivisions = subdivisions)$value
  }
}

base_integration <- function() {
  for (kernel in kernels) {
    integrate(f = kernel$fun,
              lower = kernel$support[1],
              upper = kernel$support[2])
  }
}

microbenchmark::microbenchmark(primitive_integration(), base_integration())
