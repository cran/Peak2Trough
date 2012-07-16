\name{plot.EMnlm}
\alias{plot.EMnlm}

\title{Plots seasonal variation, trend and linear predictor obtained from \code{\link{rrDynamic}} or \code{\link{rrStatic}}}

\description{Function takes as input the output from the function \code{\link{rrDynamic}} or \code{\link{rrStatic}} and plots 6 random instances of the seasonal variation and the secular trend and linear predictor.}
\usage{
\method{plot}{EMnlm}(x, ...)
}

\arguments{
  \item{x}{object of class \code{EMnlm} obtained by using \code{\link{rrDynamic}} or \code{\link{rrStatic}}.}
  \item{...}{ignored}
}

\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrDynamic}}, \code{\link{rrStatic}}
}

