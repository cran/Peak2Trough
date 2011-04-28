\name{model4}
\alias{model4}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation model 4: Seasonal variation and secular trend on log scale
}
\description{
This function simulates daily (default) or monthly counts according to a mean value structure on log scale including a seasonal variation and a secular trend.
}
\usage{
model4(RR = 1, years = 1, num.dta = 10, daily = TRUE,
       sample.size = 500, Xs, Bs, Xt, Bt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RR}{a single value. Must be greater or equal to 1}
  \item{years}{number of years to be simulated}
  \item{num.dta}{integer. Number of data sets to be simulated}
  \item{daily}{logical; if TRUE (default) daily counts are simulated, if FALSE monthly counts are simulated }
  \item{sample.size}{integer. Average number of counts per year}
  \item{Xs}{matrix. Design matrix for a seasonal variation on log scale}
  \item{Bs}{vector. Coefficient vector for a seasonal variation on log scale}
  \item{Xt}{matrix. Design matrix for a secular trend on log scale}
  \item{Bt}{vector. Coefficient vector for a secular trend on log scale}

}
\details{
Conditional on \eqn{N}, daily (default) counts, \eqn{(n_1,\ldots,n_K)}, are simulated according to a multinomial distribution with porbability parameters given by
\deqn{p_t \propto \exp(\alpha Xs\cdot Bs+Xt\cdot Bt),}{p_t propto exp(\alpha*Xs*Bs+Xt*Bt),}
\eqn{t=1,\ldots,K}, where
\deqn{\alpha=\frac{\log(RR)}{\max(Xs\cdot Bs)-\min(Xs\cdot Bs)}.}{\alpha=log(RR)/[max(Xs*Bs)-min(Xs*Bs)].}  
}
\value{           
  \item{sim.data}{list with \env{num.dta} matrices containing simulated counts}
  \item{sim.agr.data}{list with \env{num.dta} matrices containing aggregated simulated counts in case of \env{years}\eqn{>1}}. Counts are aggregated by month and year}
\references{
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
\code{\link{rr1}}, \code{\link{rr2}}, \code{\link{rr3}}
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

model4(RR=1.2,years=2,Xs=Xs,Bs=Bs,Xt=Xt,Bt=Bt)
}
\keyword{ Multimodal seasonal variation }
\keyword{ Secular trend }
\keyword{ Monte Carlo simulation }
