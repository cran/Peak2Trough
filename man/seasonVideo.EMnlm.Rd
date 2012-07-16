\name{seasonVideo.EMnlm}
\alias{seasonVideo.EMnlm}

\title{Creates video of seasonal variation component from \code{rrDynamic}}

\description{Function takes as input the output from the function \code{\link{rrDynamic}} or \code{\link{rrStatic}} and creates a video using the function \link[animation]{saveVideo}.}
\usage{
\method{seasonVideo}{EMnlm}(ss, ylim, ...)
}

\arguments{
  \item{ss}{object of class \code{EMnlm} obtained by using \code{\link{rrDynamic}} or \code{\link{rrStatic}}.}
  \item{ylim}{limits for y-axis.}
  \item{...}{inputs to \link[animation]{saveVideo}.}
}

\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrDynamic}}, \code{\link{rrStatic}}
}

