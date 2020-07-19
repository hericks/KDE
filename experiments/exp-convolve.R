grid <- 5001
sample <- 0
kernel_h_ap <- kernel_transform(rectangular,0,0.5)
ker_h <- kernel_transform(rectangular,sample,0.5)

a <- -100
b <- 100
y <- c(seq(from=sample, to= a, length.out=(as.integer(grid/2))),
       sample,
       seq(from=sample, to= b, length.out=as.integer(grid/2)))

vec <- convolve(kernel_h_ap$fun(y), rev(ker_h$fun(y)), type='open')
vec1 <- vec[1:grid]
vec2 <- vec[grid:length(vec)]


vec[as.integer(grid/2)+1]
plot(y, vec1, lwd=0.1 ,lty=1, col="red", xlim=c(-10,10))
lines(y, gaussian_function(y), lty=1)
lines(y, vec2)
