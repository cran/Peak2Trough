\name{rr2}
\alias{rr2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimating peak-to-trough ratio of a seasonal variation component modelled by Poisson regression (exact)}
\description{
Estimates the peak-to-trough ratio of a seasonal component, which may be interpreted as a relative risk with the trough as a reference group. The estimator is derived from a Poisson regression model with a single (seasonal variation) component in the linear predictor.
}
\usage{
rr2(counts, period=365, offz=NULL, plot=TRUE, origin="1960-01-01",
hemisphere=c("Northern","Southern"), x.axis=c("Months","Seasons"),
y.range=list(NULL,NULL))
}

\arguments{
  \item{counts}{a vector of observered counts. This vector may contain daily or monthly observations from at least one year.}
  \item{period}{integer stating the resolution of the year. Default is 365 indicating daily observations. Alternatives may be monthly (period=12) or weekly (period=529 observations.}
  \item{offz}{a vector of offsets for the glm call. Default is NULL.}
  \item{plot}{logical; if TRUE (default) a plot of the observed counts is made with the seasonal variation curve superimposed as well as the peak time.}
  \item{origin}{start date of study period.}
  \item{hemisphere}{character indicating whether data origins from the
  Northern (default) or Southern hemisphere. In case of Southern and
  plot=TRUE the seasonal variation curve is shifted six months.}
  \item{x.axis}{character indicating if the x axis should be divided
  into months or seasons. In case of x.axis="Season" the x axis starts
  with winter, defined as December to February in the Northern
  hemisphere and as June to August in the Sourthern hemisphere.}
  \item{y.range}{a list of two vectors giving the range of the y axes in
  the two plots.}
}
\details{
The observed counts, \eqn{n_1,\ldots,n_{K}}{n_1,\ldots,n_K}, are modelled according to the model given by

\deqn{n_t \sim \textrm{Poisson}(\lambda_t).}{n_t ~ Poisson(\lambda_t).} Using log-link we have
\deqn{\log(\lambda_t)=S_t ,}{log(\lambda_t)=S_t ,} where the seasonal component, \eqn{S_t} is given by
\deqn{S_t= \alpha\cos(\theta_t+\phi) = \gamma_1\cos(\theta_t) + \gamma_2\sin(\theta_t),}{S_t= \alphacos(\theta_t+\phi) = \gamma_1cos(\theta_t) + \gamma_2sin(\theta_t),}
for \eqn{t=1,\ldots,K}{t=1,\ldots,K}, where \eqn{\alpha \geq 0}{\alpha >= 0} is the amplitude, \eqn{\theta_t=\frac{2\pi t}{365}}{\theta_t=2\pi t/365} (if daily observations otherwise divided by 12 instead of 365) and \eqn{\phi = \arctan(\frac{\gamma_2}{\gamma_1})}{\phi = arctan(\gamma_2/\gamma_1)}.

The amplitude, \eqn{\alpha}, is estimated by

\deqn{\hat{\alpha} = \sqrt{ \gamma_1^2 + \gamma_2^2 },}{\hat{\alpha} = \sqrt \gamma_1^2 + \gamma_2^2 ,}
and the relative risk, \eqn{RR}, estimated by the peak-to-trough ratio, is then estimated by
\deqn{\hat{RR}=\exp{ 2 \hat{\alpha} }.}{RR=exp( 2 \alpha ).}

}

\value{
  \item{Model}{Formulated Poisson regression model}
  \item{RelativeRisk}{The estimated peak-to-trough ratio}
  \item{RiskDif}{The estimated peak-to-trough difference}
  \item{PeakTime}{The estimated peak time}
}
\references{
J. A. Nelder and R. W. M. Wedderburn. Generalized linear models. Journal of Royal Statistical Society 135 (1972) 370-384.

A. L. Christensen, S. Lundbye-Christensen and C. Dethlefsen. Poisson regression models outperform the geometrical model in estimating the peak-to-trough ratio of seasonal variation: a simulation study. Submitted, 2011.
}
\author{
Anette Luther Christensen
}
\note{
The code is developed in relation to and used in the manuscript:

A. L. Christensen, S. Lundbye-Christensen and C. Dethlefsen. Poisson regression models outperform the geometrical model in estimating the peak-to-trough ratio of seasonal variation: a simulation study. Submitted, 2011.
}


\seealso{
\code{\link{rr1}}, \code{\link{rr3}}, \code{\link{model1}}
}
\examples{
counts <- model1(RR=1.2,years=1,num.dta=1,sample.size=50)$sim.data[[1]]

rr2(counts)
}

\keyword{ Seasonal variation }
\keyword{ Poisson regression }
\keyword{ Peak-to-trough ratio }
