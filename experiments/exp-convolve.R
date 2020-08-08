grid <- 501
sample <- 1
kernel_h_ap <- kernel_transform(rectangular, 0, 0.5)
ker_h <- kernel_transform(rectangular, sample, 0.5)



a <- -1
b <- 1
s <- (b - a) / grid

x <- c(rev(seq(from=sample - s, to= a, length.out=(as.integer(grid/2)))),
       sample,
       seq(from=sample + s, to= b, length.out=as.integer(grid/2)))

x
vec <- s * convolve(kernel_h_ap$fun(x), rev(ker_h$fun(x)), type='open')
x_out <- seq(2*a+s, 2*b-s, len = 2*grid-1)
vec1 <- vec[1:grid]
vec2 <- vec[grid:length(vec)]

vec

vec[as.integer(grid)]
plot(x_out, vec, lwd=0.1 ,lty=1, col="red", xlim=c(-10,10))
#lines(x, triangular$fun(y))
#lines(x_out, rectangular$fun(y), lty=1)
#lines(x, vec2)
