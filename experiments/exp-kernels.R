x <- seq(from = 0, to = 10, length.out = 1000)
y <- 5
h <- 2

pal <- rainbow(11)
plot(x, rectangular((x-y)/h),
     xlim=c(0,10), ylim=c(0,1.2),
     main="Kernels", xlab="", ylab="",
     col=pal[1], type="l")
lines(x,triangular((x-y)/h), col=pal[2])
lines(x, epanechnikov((x-y)/h), col=pal[3])
lines(x, biweight((x-y)/h), col=pal[4])
lines(x, triweight((x-y)/h), col=pal[5])
lines(x, tricube((x-y)/h), col=pal[6])
lines(x, gaussian((x-y)/h), col=pal[7])
lines(x, cosine((x-y)/h), col=pal[8])
lines(x, logistic((x-y)/h), col=pal[9])
lines(x, sigmoid_fu((x-y)/h), col=pal[10])
lines(x, silverman((x-y)/h), col=pal[11])
leg_txt <- c("rectangular", "triangular", "epanechnikov",
             "biweight","triweight", "tricube", "gaussian", "cosine", "logistic",
             "sigmoid_fu", "silverman")
legend("topright", leg_txt, fill=pal[1:11])

