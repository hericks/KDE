x <- seq(from = -3, to = 8, length.out = 1000)
y <- 1:5
h <- 1

kernels <- list(expr(rectangular), expr(triangular), expr(epanechnikov),
                expr(biweight), expr(triweight), expr(tricube),
                expr(gaussian), expr(cosine), expr(logistic),
                expr(sigmoid_fu), expr(silverman))
kernel_density_est <- lapply(kernels, function(k) kde(y, h, k))

pal <- rainbow(11)
plot(x, kernel_density_est[[1]](x),
     xlim=c(-3,8), ylim=c(0,1),
     main="Kernels", xlab="", ylab="",
     col=pal[1], type="l")
lines(x,rep(0, length(x)))

for(i in 2:length(kernel_density_est)){
  lines(x, kernel_density_est[[i]](x), col=pal[i])
}
leg_txt <- c("rectangular", "triangular", "epanechnikov",
             "biweight","triweight", "tricube", "gaussian", "cosine", "logistic",
             "sigmoid_fu", "silverman")
legend("topright", leg_txt, fill=pal[1:11])
