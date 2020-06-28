kernelTransform <- function(kernel, x0, h) {
  # TODO: Check conditions on arguments.
  force(kernel)
  force(x0)
  force(h)

  function(x) {
    kernel((x-x0)/h)/h
  }
}
