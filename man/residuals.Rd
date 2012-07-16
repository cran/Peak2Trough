\name{residuals.EMnlm}
\alias{residuals.EMnlm}

\title{Residuals of a state space model}

\description{Calculation of residuals of a state space model as outlined by Jorgensen et al.}
\usage{
\method{residuals}{EMnlm}(ss)
}

\arguments{
  \item{ss}{an object os class \code{EMnlm} obtained from \code{\link{rrDynamic}} or \code{\link{rrStatic}}.}
}
\details{
The observed counts, \eqn{n_1,\ldots,n_{K}}{n_1,\ldots,n_K}, are modelled according to the dynamic generalised linear model given by
    \deqn{y_t|\lambda_t     \sim \textrm{Poisson}(\textrm{offset}_t \lambda_t)}{y_t|\lambda_t     ~ Poisson(offset_t \lambda_t)}
    \deqn{\eta_t            = \log(\lambda_t) = F_t^T  \theta_t}{\eta_t          = log(\lambda_t) = F_t^T  \theta_t}
    \deqn{\theta_t         = G_t   \theta_{t-1} + w_t, w_t \sim N(0,W)}{\theta_t = G_t  * \theta_{t-1} + w_t, w_t ~ N(0,W)}
    for \eqn{t=1,\ldots,K}{t=1,...,K}. The matrices \eqn{F_t}, \eqn{G_t},
    and \eqn{W_t} may depend on a parameter vector \eqn{\phi}.
    The initialization is given as
\deqn{\theta_0 \sim N(m_0,C_0).}{\theta_0 ~ N(m_0,C_0).}
}
\value{
Object of class \code{residuals} containing five objects:

  \item{FyRes}{filter residuals of the observation model}
  \item{SyRes}{smoother residuals of the observation model}
  \item{FthetaRes}{filter residuals of the latent model}
  \item{SthetaRes}{smoother residuals of the latent model}
  \item{misc}{object containing:}
  \item{misc$ss}{as input}
  \item{misc$FITy}{fitted (filtered) values of the observation model}
  \item{misc$FITtheta}{fitted (filtered) values of the latent model}
  \item{misc$FITsy}{fitted (smoothed) values of the observation model}
  \item{misc$FITstheta}{fitted (smoothed) values of the latent model}
}
\references{
Jorgensen, Lundbye-Christensen, Song and Sun. A state space model for multivariate longitudinal count data. Biometrika. 86 (1999) 169-181
}
\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrStatic}}, \code{\link{rrDynamic}}
}

\keyword{ Seasonal variation }
\keyword{ State Space Model }
\keyword{ Peak-to-trough ratio }

