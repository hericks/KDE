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
stepsize = 0.01
max_length = Inf
fixed = 1000

# Error Analysis
primitive_abs_errors <- double(length(kernels))
for (i in seq_along(kernels)) {
  kernel <- kernels[[i]]
  I <- integrate_primitive(integrand = kernel$fun,
                           lower = kernel$support[1],
                           upper = kernel$support[2],
                           stepsize = stepsize,
                           max_length = max_length,
                           fixed = fixed)

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
                        stepsize = stepsize,
                        max_length = max_length,
                        fixed = fixed)
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
