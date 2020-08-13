kernel <- gaussian
subdivisions = 250L

ker_h_min <- kernel_transform(kernel, 0, 0.001, subdivisions)
ker_h <- kernel_transform(kernel, 0, 0.5, subdivisions)

# crashing integral
integrate(
  function(x) {
    (ker_h_min$fun(x) - ker_h$fun(x))^2
  },
  lower = min(ker_h_min$support[1], ker_h$support[1]),
  upper = max(ker_h_min$support[2], ker_h$support[2]),
  subdivisions = subdivisions
)

# alternative

integrate(
  function(x) {
    ker_h_min$fun(x)^2
  },
  lower = min(ker_h_min$support[1], ker_h$support[1]),
  upper = max(ker_h_min$support[2], ker_h$support[2]),
  subdivisions = subdivisions
)

integrate(
  function(x) {
    (ker_h_min$fun(x) * ker_h$fun(x))
  },
  lower = min(ker_h_min$support[1], ker_h$support[1]),
  upper = max(ker_h_min$support[2], ker_h$support[2]),
  subdivisions = subdivisions
)

integrate(
  function(x) {
    ker_h$fun(x)^2
  },
  lower = min(ker_h_min$support[1], ker_h$support[1]),
  upper = max(ker_h_min$support[2], ker_h$support[2]),
  subdivisions = subdivisions
)
