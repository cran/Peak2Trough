est.phi <- function(ss, offset, dynamic=TRUE, ...){
options(warn=-1)

# Converting model from sspir notation to KFAS notation using internal function
tvar <- c(ss$n, 1, 1, 1)
Zt <- sspir2kfas(ss$Fmat, dim=c(ss$d, ss$p, tvar[1]), x=ss$x, phi=ss$phi)
Tt <- sspir2kfas(ss$Gmat, dim=c(ss$p, ss$p, tvar[2]), x=ss$x, phi=ss$phi)
Qt <- sspir2kfas(ss$Wmat, dim=c(ss$p, ss$p, tvar[4]), x=ss$x, phi=ss$phi)

# Changing initial values to fit to KFAS
a1 <- Tt[,,1] %*% matrix(ss$m0, ncol=1)
P1 <- (Tt[,,1] %*% ss$C0 %*% t(Tt[,,1]) + Qt[,,1])

# specify DGLM in KFAS
m3kfas <- KFAS::SSModel(y=ss$y, Z=Zt, T=Tt, R=diag(ss$p), Q=Qt, a1=a1, P1=P1, u=offset, distribution="Poisson")

# Log likelihood function to be minimized (maximized)
likf <- function(par, nsim, ...){
phi <- x2phi(par, ...) #transforming hyper parametervector to be defined on R^17
#updating covariance matrix of the state model 
try(# necessary in case of invalid hyper parametervectors -> numerical instability in nlm()
m3kfas$Q[,,1] <- as.matrix(nearPD(phi2W(phi, ...) )$mat)
#.sidste <<- -logLik(m3kfas, nsim=nsim), silent = TRUE)
#-logLik(m3kfas, nsim=nsim)
, silent=TRUE)
if(class(.Last.value)=="try-error"){ return(m3kfas$loglik) } #when invalid hyper parametervectors log likelihood value assigned to last valid value, hence totally flat likelihood function (consequences not investigated)
m3kfas$loglik <- -logLik(m3kfas, nsim=nsim)
#assign(".sidste", -logLik(m3kfas, nsim=nsim), envir=.GlobalEnv)
return(m3kfas$loglik)
}

likfstat <- function(par, nsim, ...){
par <- c(par,rep(0,2*ss$S))
phi <- x2phi(par, ...) #transforming hyper parametervector to be defined on R^17
#updating covariance matrix of the state model 
try(# necessary in case of invalid hyper parametervectors -> numerical instability in nlm()
m3kfas$Q[,,1] <- as.matrix(nearPD(phi2W(phi, ...) )$mat), silent=TRUE)
if(class(.Last.value)=="try-error"){ return(m3kfas$loglik) } 
#.sidste <<- -logLik(m3kfas, nsim=nsim), silent = TRUE)
#assign(".sidste", -logLik(m3kfas, nsim=nsim), envir=.GlobalEnv)
#-logLik(m3kfas, nsim=nsim)
m3kfas$loglik <- -logLik(m3kfas, nsim=nsim)
#when invalid hyper parametervectors log likelihood value assigned to last valid value, hence totally flat likelihood function (consequences not investigated)
return(m3kfas$loglik)
}

# obtaining estimate of hyper parametervector using nlm()
if(!dynamic){
start <- phi2x(ss$phi, ...)[1:(ss$T*(ss$T+1)/2)]
time.nlm <- system.time(fit1 <- nlm(p=start, f=likfstat, nsim=200, gradtol=1e-12, print.level=0,  ...))
fit1$estimate <- c(fit1$estimate, rep(0,2*ss$S))
}else{
start <- phi2x(ss$phi, ...)
time.nlm <- system.time(fit1 <- nlm(p=start, f=likf, nsim=200, gradtol=1e-12, print.level=0, ...))
}

m4kfas <- m3kfas
# udpdating Q
m4kfas$Q[,,1] <- as.matrix(nearPD(phi2W(x2phi(fit1$estimate, ...), ...))$mat)
# calculates approximating Gaussian model
out <- KFAS::approxSSM(m4kfas)

# creates Gaussian model to be smoothed
m5kfas <- m4kfas
m5kfas$distribution <- "Gaussian"
m5kfas$y <- ss$ytilde <- out$y
m5kfas$H <- ss$vtilde <- out$H

outg <- KFAS::KFS(m5kfas, smoothing="state")

# Backtracking filtered means, m, and covariances, C
      ss$mf <- t( (solve(Tt[,,1:tvar[2]]) %*% outg$a)[,-c(1)] )
      ss$Cf <- list()      
      for(j in 1:ss$n){
          ss$Cf[[j]] <- solve(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]) %*% (outg$P[,,j+1] - Qt[,,(tvar[4]!=1)*(j+2)+(tvar[4]==1)]) %*% solve(t(Tt[,,(tvar[2]!=1)*(j+2)+(tvar[2]==1)]))
      }

      # Backtracking smoothed means, m.tilde, and covariances, C.tilde
      ss$m <- t(outg$alphahat)
      ss$C <- list()
      for(j in 1:ss$n){
          ss$C[[j]] <- outg$V[,,j]
      }
ss$phi <- x2phi(fit1$estimate, ...)
ss$Wmat <- function(tt, x, phi) phi2W(ss$phi, ...)
ss$loglik <- -fit1$minimum
options(warn=0)
list(ss=ss, nlm.fit=fit1)
}

