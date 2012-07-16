\name{plot.residuals}
\alias{plot.residuals}

\title{Plots residuals from \code{\link{residuals.EMnlm}}}

\description{Function takes as input the output from the function \code{\link{residuals.EMnlm}} and plots all residuals.}
\usage{
\method{plot}{residuals}(
	resid, 
	path  = NULL, 
	mtextx=c(expression(T[t-1]), expression(a[t-1]),
		expression(alpha[list(1,t-1)]), expression(beta[list(1,t-1)]),
		expression(alpha[list(2,t-1)]), expression(beta[list(2,t-1)]),
		expression(tilde(T)[t-1]), expression(tilde(a)[t-1]),
		expression(tilde(alpha)[list(1,t-1)]), expression(tilde(beta)[list(1,t-1)]),
		expression(tilde(alpha)[list(2,t-1)]), expression(tilde(beta)[list(2,t-1)])), 
	mtexty=c(expression(T[t]), expression(a[t]), expression(alpha[list(1,t)]),
		expression(beta[list(1,t)]), expression(alpha[list(2,t)]),
		expression(beta[list(2,t)]), expression(tilde(T)[t]),
		expression(tilde(a)[t]), expression(tilde(alpha)[list(1,t)]),
		expression(tilde(beta)[list(1,t)]), expression(tilde(alpha)[list(2,t)]),
		expression(tilde(beta)[list(2,t)])), 
	burn.in=1, ...)
}

\arguments{
  \item{resid}{object of class \code{residuals} obtained by using \code{\link{residuals.EMnlm}} or \code{\link{rrStatic}}.}
  \item{path}{path where to save plots.}
  \item{mtextx}{vector with text for the x axes using the \code{\link[base]{expression}} function.}
  \item{mtexty}{vector with text for the y axes using the \code{\link[base]{expression}} function.}
  \item{burn.in}{number starting residuals to ignore.}
  \item{...}{options to \code{\link[grDevices]{pdf}} function.}
}

\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrDynamic}}, \code{\link{rrStatic}}, \code{\link{residuals.EMnlm}}
}

