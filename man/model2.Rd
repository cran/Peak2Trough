\name{model2}
\alias{model2}
\title{
Simulation model 2: Unimodal sinusoid on original scale including a secular trend
}
\description{
This function simulates daily (default) or monthly counts according to a mean value structure on original scale including a unimodal sinusoid and a trend on log scale.
}
\usage{
model2(RR = 1, years = 1, num.dta = 10, daily = TRUE, 
       sample.size = 500, Xt, Bt)
}
\arguments{
  \item{RR}{a single value. Must be greater or equal to 1}
  \item{years}{number of years to be simulated}
  \item{num.dta}{integer. Number of data sets to be simulated}
  \item{daily}{logical; if TRUE (default) daily counts are simulated, if FALSE monthly counts are simulated }
  \item{sample.size}{integer. Average number of counts per year}
  \item{Xt}{matrix. Design matrix for a secular trend on log scale}
  \item{Bt}{vector. Coefficient vector for a secular trend on log scale}
}
\details{
Conditional on \eqn{N}, daily (default) counts, \eqn{(n_1,\ldots,n_K)}{(n_1,\ldots,n_K)}, are simulated according to a multinomial distribution with porbability parameters given by
\deqn{p_t \propto [1 + \alpha\cos(\theta_t - \theta^*)]\exp(XtBt),}{p_t propto [1 + \alpha cos(\theta_t - \theta*)]exp(XtBt),}
\eqn{t=1,\ldots,K}, where
\deqn{\alpha=\frac{RR-1}{RR+1},}{\alpha=(RR-1)/(RR+1),}  and \eqn{\theta^* \in [0,2\pi]}{0 <= \theta* <= 2\pi} is chosen randomly. 
}
\value{           
  \item{sim.data}{list with \env{num.dta} matrices containing simulated counts}
  \item{sim.agr.data}{list with \env{num.dta} matrices containing aggregated simulated counts in case of \env{years}\eqn{>1}. Counts are aggregated by month and year}
}
\references{
J. H. Edwards. The recognition and estimation of cyclic trends. Annals of Human Genetics 25 (1961) 83-86.

M. A. Brookhart and K. J. Rothman. Simple estimators of the intensity of seasonal occurrence. BMC Medical Research Methodology 8 (2008) 67-75.

A. L. Christensen, S. Lundbye-Christensen and C. Dethlefsen. Poisson regression models outperform the geometrical model in estimating the peak-to-trough ratio of seasonal variation: a simulation study. Submitted, 2011.
}
\author{Anette Luther Christensen}
\note{
The code is developed in relation to and used in the manuscript: 

A. L. Christensen, S. Lundbye-Christensen and C. Dethlefsen. Poisson regression models outperform the geometrical model in estimating the peak-to-trough ratio of seasonal variation: a simulation study. Submitted, 2011.
}
\seealso{
\code{\link{rr1}}, \code{\link{rr2}}, \code{\link{rr3}}
}
\examples{
years <- 2
# Linear secular trend on log scale
Xt <- cbind(1,1:(365*years))
Bt <- matrix(c(20,0.005),ncol=1)
model2(RR=1.2,years=2,Xt=Xt,Bt=Bt)
}

\keyword{ Seasonal variation }
\keyword{ Secular trend }
\keyword{ Monte Carlo simulation }
