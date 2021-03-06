% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/builtin_kernels.R
\docType{data}
\name{triweight}
\alias{triweight}
\title{Triweight Kernel}
\format{
An object of S3 class \code{\link{Kernel}} with named entries
\itemize{
\item \bold{\code{fun}}: the triweight function
\item \bold{\code{support}}: \code{c(-1,1)}
\item \bold{\code{subdivisions}}: \code{1000L}.
}
}
\usage{
triweight
}
\description{
The triweight kernel is S3 object of class \code{Kernel}
provided by the \code{KDE} package.
}
\seealso{
\code{\link{Kernel}} for more information about kernels.

Other kernels: 
\code{\link{biweight}},
\code{\link{cosine}},
\code{\link{epanechnikov}},
\code{\link{gaussian}},
\code{\link{logistic}},
\code{\link{rectangular}},
\code{\link{sigmoid}},
\code{\link{silverman}},
\code{\link{triangular}},
\code{\link{tricube}}
}
\concept{kernels}
\keyword{datasets}
