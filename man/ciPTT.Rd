\name{ciPTT}
\alias{ciPTT}

\title{Confidence intervals of peak-to-trough ratio}

\description{Estimates confidence intervals for the peak-to-trough ratio of a seasonal variation component estimated by \code{\link{rr2}} or \code{\link{rr3}} by use of block-bootstrapping.}
\usage{
ciPTT(y, nboot, b.length=10, estimator, alpha=c(0.025, 0.975))
}

\arguments{
  \item{y}{vector of counts.}
  \item{nboot}{number of bootstraps used.}
  \item{b.length}{length of block. Default 10.}
  \item{estimator}{function of estimator.}
  \item{alpha}{ends of confidence intervals. Default 0.025 and 0.975.}
}
\references{
Efron and Tibshirani. An introduction to the Bootstrap. Chapman & Hall/CRC, 1st edition, 1998.

Edwards, JH. The recognition and estimation of cyclic trends. Annals of Human Genetics. 25:83-87, 1961.
}
\author{
Anette Luther Christensen
}
\seealso{
\code{\link{rr2}}, \code{\link{rr3}}
}
\examples{
est <- function(y){ rr2(y, period=12, plot=FALSE)$RR }

## Counts from Edwards (1961) ##
y <- c(10, 19, 18, 15, 11, 13, 7, 10, 13, 23, 15, 22)

fit <- rr2(y, period=12, plot=FALSE)

ciPTT(fit$Model$y, nboot=1000, estimator=est, b.length=2)$confpoints
}
