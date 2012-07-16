\name{rrDynamic}
\alias{rrDynamic}

\title{Estimating peak-to-trough ratio of a dynamic seasonal variation component modelled by state space models}

\description{Estimates the changing peak-to-trough ratio of a dynamic seasonal variation component modelled by a state space model. The function fits a specific structured state space model, assuming data being Poisson distributed and the latent process to be Gaussian. The linear predictor is assumed being a superposition of a dynamic secular trend and a dynamic seasonal variation. The function estimates the covariance matrix of the latent process, \eqn{W}{W}, assumed to be constant, by maximising the likelihood function by use of initially the EM algorithm and switching to the \code{nlm} routine according options given by the user.}
\usage{
rrDynamic(counts, 
          risktime, 
          offset   = 1, 
          dates, 
          phi, 
          max.iter = 1000, 
          epsilon  = 0.1, 
          T        = 2, 
          S        = 2, 
          C0       = 100*diag(T+2*S), 
          m0       = matrix(rep(0, T+2*S), nrow=1), 
          period   = 52, 
          subset   = 1:length(counts),
          V        = 1,
	  estimate = TRUE,
	  ...)
}

\arguments{
  \item{counts}{a vector of observered counts.}

  \item{risktime}{a vector of risktime for each observation. Same length a \code{y}.}

  \item{offset}{constant value for all observations.}

  \item{dates}{a vector of class \code{date} specifying the date for each observation.}

  \item{phi}{a vector of starting hyper parameters.}

  \item{max.iter}{number of iterations before switching from the EM algorithm to the \code{nlm} routine for optimizing the log likelihood function.}

  \item{epsilon}{convergence criterion for the EM algorithm. If obtained before the number of iterations exceeds \code{max.iter}, the optimisation if switched to the \code{nlm} routine.}

  \item{T}{dimension of the secular trend. If \code{T=1} the trend is modelled as a random walk, if \code{T=2} the trend is modelled as a local linear trend, if \code{T=p} the trend is modelled as a \eqn{p}{p}th order local polynomial growth.}

  \item{S}{number of sinusoids in the seasonal variation: \eqn{\displaystyle\sum_{i=1}^S a_{i,t}\cos(i\theta_t)+b_{i,t}\sin(i\theta_t)}{\sum^S a_it\cos(i\theta_t)+b_it\sin(i\theta_t)}, \eqn{\theta_t=\frac{2\pi t}{\texttt{period}}}{\theta_t=2\pi t/\code{period}}.}

  \item{C0}{covariance matrix of prior distribution of the latent process.}

  \item{m0}{mean vector of prior distribution of the latent process.}

  \item{period}{the period of seasonality, eg 365 if daily observations, 52 if weekly observations, or 24 if hourly observations and the circadian variation is of interest.}

  \item{subset}{a vector if only a subset of the data should be analysed.}

  \item{V}{the variance of the observation model.}

  \item{estimate}{logical. TRUE (default) if estimation of hyperparameters should be performed. Otherwise only smoothing is performed.}
  \item{...}{options passed to \code{\link[stats]{nlm}}.}
}
\details{
The observed counts, \eqn{n_1,\ldots,n_{K}}{n_1,\ldots,n_K}, are modelled according to the dynamic generalised linear model given by

\deqn{
y_t|\lambda_t \sim \textrm{Poisson}(\textrm{offset}_t \lambda_t)
}{
y_t|\lambda_t ~ Poisson(offset_t \lambda_t)
}

\deqn{
\eta_t = \log(\lambda_t) = F_t^\top  \theta_t
}{
\eta_t = log(\lambda_t) = F_t^T  \theta_t
}

\deqn{
\theta_t = G_t \theta_{t-1} + w_t, w_t \sim N(0,W)
}{
\theta_t = G_t  * \theta_{t-1} + w_t, w_t ~ N(0,W)
}
for \eqn{t=1,\ldots,K}{t=1,...,K}. The matrices \eqn{F_t}, \eqn{G_t}, and \eqn{W_t} may depend on a parameter vector \eqn{\phi}. The initialization is given as

\deqn{
\theta_0 \sim N(m_0,C_0).
}{
\theta_0 ~ N(m_0,C_0).
}

The latent process consists of two components, a secular trend and a seasonal variation, ie 
\deqn{
\theta_t^\top = [\theta_t^{(T)\top} \; \theta_t^{(S)\top}].
}{
\theta_t^T = [T.\theta_t^T S.\theta_t^T].
}
The trend is specified as a polynomial of order \eqn{p}, denoted \eqn{q(t)}, hence
\deqn{
\theta_t^{(T)\top} = [q(t) \; q'(t) \; \cdots \; q^{(p)}(t)]
}{
T.\theta_t^T = [q(t) q'(t) ... q^p(t)]
}
and the seasonal variation is specified as a sum of \eqn{k} sinusoids, hence

\deqn{
\theta_t^{(S)\top} = [cs_{1,t} \; \ldots \; cs_{k,t}],
}{
S.\theta_t^T = [cs_{1,t} ... cs_{k,t}],
}
where \eqn{cs_{i,t}=[\cos(i\frac{2\pi}{\texttt{period}}t) \;
  \sin(i\frac{2\pi}{\texttt{period}}t)]}{cs_{i,t}=[cos(it*(2\pi/\code(period)))
  sin(it*(2\pi/\code(period)))]} for \eqn{i=1,\ldots,k}{i=1,...,k}.

The design matrix of the latent process, \eqn{G_t} becomes




}
\value{
Object of class \code{EMnlm} containing three objects:

\code{em} object is \code{NULL} if \code{estimate=FALSE}
  \item{ss}{object of class \code{\link[sspir]{SS}} with updated values according to the EM algorithm}
  \item{fit}{object containing:}
  \item{fit$time}{computation time for the EM algorithm}
  \item{fit$iterations}{number of iterations upon convergence}
  \item{fit$convergence}{logical; \code{TRUE} if convergence reached before number of iterations exceeds \code{max.iter}}
  \item{fit$loglik}{log likelihood value upon convergence}
  \item{fit$estimate}{parameter estimate upon convergence}

\code{nlm} object is \code{NULL} if \code{estimate=FALSE}
  \item{ss}{object of class \code{\link[sspir]{SS}} with updated values according to the \code{\link[stats]{nlm}} routine}
  \item{AIC}{a list of dimension two: AICs calculated as outlined by both Gersch & Kitagawa (1983) and Durbin & Koopman (2001)}
  \item{fit}{object containing:}
  \item{fit$time}{computation time for the nlm routine}
  \item{fit$iterations}{number of iterations upon convergence}
  \item{fit$convergence}{code for convergence of the \code{nlm} routine, see \code{\link[stats]{nlm}} for details}
  \item{fit$loglik}{log likelihood value upon convergence}
  \item{fit$estimate}{parameter estimate upon convergence}


\code{misc}
  \item{S}{as input}
  \item{T}{as input}
  \item{V0}{starting value of \code{V}}
  \item{phi0}{starting value of \code{phi}}
  \item{m0}{as input}
  \item{C0}{as input}
  \item{period}{as input}
  \item{subset}{as input}
  \item{data}{data frame containing \code{counts}, \code{risktime}, and \code{dates}}
  \item{offset}{as input}
  \item{epsilon}{as input}
  \item{max.iter}{as input}

\code{results}
  \item{ptt, ptt.dif, time.min, time.max}{Peak-to-trough ratio}
  \item{ptt.dif}{Peak-to-trough difference}
  \item{time.min}{time for trough}
  \item{time.max}{time for peak}
}
\references{
West and Harrison. Bayesian Forecasting and Dynamic Models. 2nd, Springer, 1997

Gersch and Kitagawa. The Prediction of Time Series With Trends and Seasonalities. Journal of Business & Economic Statistics. 1 (1983) 253-264

Durbin and Koopman. Time Series Analysis by STate Space Methods. 1st, Oxford University Press, 2001
}
\author{
Anette Luther Christensen
}

\seealso{
\code{\link{rrStatic}}
}

\keyword{ Seasonal variation }
\keyword{ State Space Model }
\keyword{ Peak-to-trough ratio }

