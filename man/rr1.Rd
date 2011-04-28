\name{rr1}
\alias{rr1}

\title{Estimating peak-to-trough ratio of a seasonal variation component modelled by the geometrical model}

\description{Estimates the peak-to-trough ratio of a seasonal component, which may be interpreted as a relative risk with the trough as a reference group. The estimator is derived by Brookhart and Rothman (2008).}
\usage{
rr1(counts, plot=TRUE)
}

\arguments{
  \item{counts}{a vector of observered (aggregated) monthly counts. This vector must contain 12 observations.}
  \item{plot}{logical; if TRUE a plot of the observed counts is made with the seasonal variation curve superimposed as well as the peak time.}
}
\details{
The observed counts, \eqn{n_1,\ldots,n_{12}}{n_1,\ldots,n_12}, are modelled according to the model given by

\deqn{n_i = \frac{N}{12}  [1+\alpha\cos(\theta_i-\theta^*)],}{n_i = N/12 * [1+\alphacos(\theta_i-\theta*)],}
for \eqn{i=1,\ldots,12}{i=1,\ldots,12}, where \eqn{N=\sum n_i}, \eqn{0 <= \alpha < 1} is the amplitude, \eqn{\theta_i=\frac{2\pi i}{12}}{\theta_i=2\pi i/12} and \eqn{\theta^*\in [0,2\pi]}{0 <= \theta* <= 2\pi}  indicates the time point with highest count. 

The amplitude, \eqn{\alpha}, is estimated by

\deqn{\hat{\alpha} = 2\sqrt{ \frac{D^212^2-Nf}{N(N-1)} },}{\alpha = 2\sqrt (D^2*12^2-N*f)/(N(N-1)) ,} 
where 
\deqn{D=\sqrt{ (12^{-1} \sum n_i\cos(\theta_i))^2 + (12^{-1} \sum n_i\sin(\theta_i))^2 }}{D=\sqrt{ (1/12 \sum n_i\cos(\theta_i))^2 + (1/12 \sum n_i\sin(\theta_i))^2 }}
and 
\deqn{f = \frac{D^212^2/N}{1+D^212^2/N},}{f = (D^2*12^2/N)/(1+D^2*12^2/N),}
the latter being a correction term to ensure valid values of the amplitude.


The relative risk, \eqn{RR}, estimated by the peak-to-trough ratio, is then estimated by

\deqn{\hat{RR} = \frac{1+\hat{\alpha}}{1-\hat{\alpha}}.}{RR = (1+\alpha)/(1-\alpha).}

}
\value{
  \item{RelativeRisk}{The estimated peak-to-trough ratio}
  \item{PeakTime}{The estimated peak time}
}
\references{
J. H. Edwards. The recognition and estimation of cyclic trends. Annals of Human Genetics 25 (1961) 83-86.

M. A. Brookhart and K. J. Rothman. Simple estimators of the intensity of seasonal occurrence. BMC Medical Research Methodology 8 (2008) 67-75.

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
\code{\link{rr2}}, \code{\link{rr3}}, \code{\link{model1}}
}
\examples{
counts <- model1(RR=1.2,years=1,daily=FALSE,num.dta=1)$sim.data[[1]]

rr1(counts)
}
\keyword{ Seasonal variation }
\keyword{ Geometrical model }
\keyword{ Peak-to-trough ratio }

