\name{ciPTT.EMnlm}
\alias{ciPTT.EMnlm}

\title{Confidence intervals of peak-to-trough ratio}

\description{Estimates confindece intervals for the peak-to-trough ratio of a seasonal variation component by use of block-bootstrapping.}
\usage{
\method{ciPTT}{EMnlm}(m, nboot, b.length=10, estimator, alpha=c(0.025, 0.975))
}

\arguments{
  \item{m}{object of class \code{EMnlm} obtained by using \code{\link{rrDynamic}} or \code{\link{rrStatic}}.}
  \item{nboot}{number of bootstraps used.}
  \item{b.length}{length of block. Default 10.}
  \item{estimator}{function of estimator.}
  \item{alpha}{ends of confidence intervals. Default 0.025 and 0.975.}
}
\references{
Efron and Tibshirani. An introduction to the Bootstrap. Chapman & Hall/CRC, 1st edition, 1998.
}
\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrDynamic}}, \code{\link{rrStatic}}
}

