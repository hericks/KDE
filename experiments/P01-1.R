test313123 <- "Test"
test <- "Test"
  test2 <- "Test"
test2 <- "Test"
test2 <- "Test"
test2 <- "Test"

1:10
lsq(matrix(letters[1:6] , nrow=2), 1:3)
lsq(matrix(1:6, nrow=3), list(1,2,3))
lsq(matrix(1:6, nrow=3), 1:4)
lsq(matrix(1:6, nrow=3), matrix(1:1, nrow=1))
lsq(matrix(1:6, nrow=3), matrix(1:1, nrow=3))
lsq(matrix(double(0), nrow=0, ncol=0), matrix(double(0), nrow=0, ncol=0))
lsq(matrix(1:6, nrow=3), c(1,NA,3))

lsq <- function(X, y) {
  # TODO
  A <- t(X) %*% X
  solve(A, t(X) %*% y)
}

lsq(matrix(1:6, nrow=3), 1:3)
lsq(matrix(runif(6), nrow=3), matrix(runif(3), ncol=1))
lsq(matrix(letters[1:6] , nrow=2), 1:3)


lsq(matrix(1:6, nrow=2), matrix(1:1, nrow=1))
lsq(matrix(1:6, nrow=3), matrix(1:1, nrow=3))
lsq(matrix(double(0), nrow=0, ncol=0), matrix(double(0), nrow=0, ncol=0))
lsq(matrix(1:6, nrow=3), c(1,NA,3))
lsq(matrix(c(1:5, NA), nrow=3), 1:3)
lsq(matrix(c(1:5, NA), ncol=3), 1:10)
lsq(matrix(c(1,1,2,1,1,2), nrow=3), 1:3)

lsq <- function(X, y) {
  # TODO
  A <- t(X) %*% X
  solve(A, t(X) %*% y)
}
