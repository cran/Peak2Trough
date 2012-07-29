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
, silent=TRUE)
if(class(.Last.value)=="try-error"){ return(m3kfas$loglik) } 
m3kfas$loglik <- -logLik(m3kfas, nsim=nsim)
return(m3kfas$loglik)
}

likfstat <- function(par, nsim, ...){
par <- c(par, length(ss$phi)-(ss$T*(ss$T+1)/2))
phi <- x2phi(par, ...) #transforming hyper parametervector to be defined on R^17
#updating covariance matrix of the state model 
try(# necessary in case of invalid hyper parametervectors -> numerical instability in nlm()
m3kfas$Q[1:ss$T,1:ss$T,1] <- as.matrix(nearPD(phi2W(phi, ...)[1:ss$T,1:ss$T] )$mat), silent=TRUE)
if(class(.Last.value)=="try-error"){ return(m3kfas$loglik) } 
m3kfas$loglik <- -logLik(m3kfas, nsim=nsim)
return(m3kfas$loglik)
}

# obtaining estimate of hyper parametervector using nlm()
if(!dynamic){
start <- phi2x(ss$phi, ...)[1:(ss$T*(ss$T+1)/2)]
time.nlm <- system.time(fit1 <- nlm(p=start, f=likfstat, typsize=abs(start), nsim=200, steptol=min(abs(start)), stepmax=max(abs(start)), ...))
fit1$estimate <- c(fit1$estimate, rep(0,length(ss$phi)-(ss$T*(ss$T+1)/2)))
m4kfas <- m3kfas
# udpdating Q
m4kfas$Q[1:ss$T, 1:ss$T, 1] <- phi2W(x2phi(fit1$estimate, ...), ...)[1:ss$T, 1:ss$T]
}else{
start <- phi2x(ss$phi, ...)
time.nlm <- system.time(fit1 <- nlm(p=start, f=likf, nsim=200, typsize=abs(start), steptol=min(abs(start)), stepmax=max(abs(start)), ...))
m4kfas <- m3kfas
# udpdating Q
m4kfas$Q[,,1] <- phi2W(x2phi(fit1$estimate, ...), ...)
}

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
if(dynamic){ss$phi <- x2phi(fit1$estimate, ...)}else{ss$phi <- c(x2phi(fit1$estimate, ...)[1:(ss$T*(ss$T+1)/2)], rep(0,length(ss$phi)-(ss$T*(ss$T+1)/2)))}
ss$Wmat <- function(tt, x, phi) phi2W(ss$phi, ...)
ss$loglik <- logLik(m5kfas, nsim=200)
options(warn=0)
list(ss=ss, nlm.fit=fit1)
}

