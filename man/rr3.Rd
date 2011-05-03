\name{rr3}
\alias{rr3}

\title{
Estimating peak-to-trough ratio of a seasonal variation component modelled by Poisson regression
}
\description{
Estimates the peak-to-trough ratio of a seasonal component, which may be interpreted as a relative risk with the trough as a reference group. The estimator is derived from a Poisson regression model with a seasonal variation component and a secular trend in the linear predictor.
}
\usage{
rr3(counts, period=365, offz=NULL, plot=TRUE, origin="1960-01-01",
hemisphere=c("Northern","Southern"), x.axis=c("Months","Seasons"),
y.range=list(NULL,NULL,NULL))
}

\arguments{
  \item{counts}{a vector of observered counts. This vector may contain daily or monthly observations from at least one year.}
  \item{period}{integer stating the resolution of the year. Default is 365 indicating daily observations. Alternatives may be monthly (period=12) or weekly (period=529 observations.}
  \item{offz}{a vector of offsets for the glm call. Default is NULL.}
  \item{plot}{logical; if TRUE (default) a three plots are made. The first of the observed counts with the linear predictor superimposed. The second of the exponentiated secular trend and the third of the detrended seasonal variation as well as the peak time.}
  \item{origin}{start date of study period.}
  \item{hemisphere}{character indicating whether data origins from the
  Northern (default) or Southern hemisphere. In case of Southern and
  plot=TRUE the seasonal variation curve is shifted six months.}
  \item{x.axis}{character indicating if the x axis should be divided
  into months or seasons. In case of x.axis="Season" the x axis starts
  with winter, defined as December to February in the Northern
  hemisphere and as June to August in the Sourthern hemisphere.}
  \item{y.range}{a list of three vectors giving the range of the y axes in
  the three plots.}
}
\details{
The observed counts, \eqn{n_1,\ldots,n_{K}}{n_1,\ldots,n_K}, are modelled according to the model given by

\deqn{n_t \sim \textrm{Poisson}(\lambda_t).}{n_t ~ Poisson(\lambda_t).} Using log-link we model a seasonal variation and a secular trend component
\deqn{\log(\lambda_t)=S_t + T_t .}{log(\lambda_t)=S_t + T_t .} The seasonal component, \eqn{S_t}, is a superposition of four sinusoids of different frequecies each given by
\deqn{S_{jt}= \alpha_j\cos(j\theta_t+\phi_j) = \gamma_{j,1}\cos(j\theta_t) + \gamma_{j,2}\sin(j\theta_t),}{S_{jt}= \alpha_jcos(j\theta_t+\phi_j) = \gamma_{j,1}cos(j\theta_t) + \gamma_{j,2}sin(j\theta_t),}
for \eqn{j=1,\ldots,4}{j=1,\ldots,4}, where \eqn{\alpha_j \geq 0}{\alpha_j >= 0} is the amplitude of the \eqn{j}th sinusoid, \eqn{\theta_t=\frac{2\pi t}{365}}{\theta_t=2\pi t/365} (if daily observations otherwise divided by 12 instead of 365) and \eqn{\phi_j = \arctan(\frac{\gamma_{j,2}}{\gamma_{j,1}})}{\phi_j = arctan(\gamma_{j,2}/\gamma_{j,1})}. The secular trend is modelled as a restricted cubic spline with five knots.

The relative risk, \eqn{RR}, is estimated by the peak-to-trough ratio, given by
\deqn{\hat{RR}=\frac{\max(\exp{ S_t })}{\min(\exp{ S_t })}.}{RR=max[exp( S_t )]/min[exp( S_t )].}

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
\code{\link{rr1}}, \code{\link{rr2}}, \code{\link{model4}}
}
\examples{
years <-2
dato <- 1:(365*years)
# Seasonal variation as a superposition of four sinsoids
# with frequencies 1 to 4 on log scale
Xs <- cbind(cos(dato*2*pi/365) , sin(dato*2*pi/365) ,
            cos(dato*4*pi/365) , sin(dato*4*pi/365) ,
            cos(dato*6*pi/365) , sin(dato*6*pi/365) ,
            cos(dato*8*pi/365) , sin(dato*8*pi/365) )
Bs <- matrix(c( 0.07  ,  0.01  ,
               -0.01  ,  0.02  ,
                0.01  ,  0.02  ,
                0.002 ,  0.003),
                ncol=1)
# Linear secular trend on log scale
Xt <- cbind(1,1:(365*years))
Bt <- matrix(c(20,0.005),ncol=1)

counts <- model4(RR=1.2,years=2,Xs=Xs,Bs=Bs,Xt=Xt,Bt=Bt,num.dta=1,sample.size=5000)$sim.data[[1]]

rr3(counts)
}
\keyword{ Seasonal variation }
\keyword{ Poisson regression }
\keyword{ Peak-to-trough ratio }

