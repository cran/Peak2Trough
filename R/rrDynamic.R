rrDynamic <- function(counts, risktime, offset=1, dates, phi, max.iter=1000, epsilon=0.1, T=2, S=2, C0=100*diag(T+2*S), m0=matrix(rep(0,T+2*S),nrow=1), period=52, subset=1:length(counts), V=1, estimate=TRUE, ...){
	
	cat("Initialising ... \n")

	n    <- length(counts)

	tvar <- c(n, 1, 1, 1)
	# creating m1:
	m1   <- ssm(counts ~ -1 +
          tvar(polytime(1:n, T-1)) +
          tvar(polytrig(1:n, period, S)),
          family = poisson,
          fit    = FALSE,
	  C0     = C0,
          m0     = m0,
          phi    = phi,
          Vmat   = function(tt,x,phi){ matrix(V, nrow=1, ncol=1) }
          )
	Wmat(m1)   <- function(tt,x,phi){ phi2W(phi, T, S) }

	if(estimate){
	# H transform m1 -> m2
	H          <- H.trans(T,S)
	m2         <- m1
	C0(m2)     <- H%*%C0%*%t(H)
	m0(m2)     <- matrix((H%*%t(m0)), nrow=1)
	phi(m2)    <- phi
	Vmat(m2)   <- function(tt, x, phi){ matrix(V,nrow=1,ncol=1) }
	Wmat(m2)   <- function(tt, x, phi){ H%*%phi2W(m2$ss$phi, T, S)%*%t(H) }
	Fmat(m2)   <- function(tt, x, phi){ t(x$x[tt,]%*%solve(H)) }
	m2$ss$S    <- S
	m2$ss$T    <- T
	
	cat("Applies the EM algorithm ... \n")

	# Kører EM på m2 -> m3
	time.em    <- system.time(
		em.out2 <- EMalgo.kfas(m2$ss, offset = risktime/offset, tvar = tvar, Wstruc = W.strc, maxiter = max.iter, epsilon = epsilon, trace = FALSE)
	)

	m3         <- m1
	phi(m3)    <- W2phi(solve(H)%*%em.out2$Wmat.est%*%t(solve(H)), T, S)
	Wmat(m3)   <- function(tt, x, phi){ phi2W(phi, T, S) }
	m0(m3)     <- matrix(solve(H)%*%t(smooth.bottom(em.out2$ss)$m0), nrow=1)
	C0(m3)     <- solve(H)%*%(smooth.bottom(em.out2$ss)$C0)%*%solve(t(H))

	em.model   <- list(ss=m3$ss, fit=list(time=time.em, iterations=em.out2$iterations, convergence=em.out2$convergence, loglik=tail(em.out2$loglik, n=1), estimate=m3$phi))

	cat("Applies the stat::nlm routine ... \n")

	# kører nlm() på m3 -> m4
	time.nlm   <- system.time(
		m4 <- est.phi(m3$ss, offset=risktime/offset, S=S, T=T, ...)
	)

	
	
	nlm.model  <- list(ss=m4$ss, fit=list(time=time.nlm, iterations=m4$nlm.fit$iterations, convergence=m4$nlm.fit$code, loglik=-m4$nlm.fit$minimum, estimate=x2phi(m4$nlm.fit$estimate, T=T, S=S)), AIC=list("K&G"=-2*(-m4$nlm.fit$minimum)+2*(T+2*S+length(phi)),"D&K"=(-2*(-m4$nlm.fit$minimum)+2*length(phi))/n))

	results    <- numeric()
	for(i in 1:nlm.model$ss$n){
		results[i] <- calc.season(nlm.model$ss, time=i, period=period, idx=T+(1:(S*2)))[1]
	}

	}else{

	m4 <- Fkfs(ss = m1$ss, tvar, offset = risktime/offset)

	results    <- numeric()
		for(i in 1:n){
			results[i] <- calc.season(m4$ss, time=i, period=period, idx=T+(1:(S*2)))[1]
		}
	em.model <- NULL
	nlm.model <- m4
	}
	cat("Finalising ... \n")
	misc       <- list(S=S, T=T, V0=V, phi0=phi, m0=m0, C0=C0, period=period, subset=subset, data=data.frame(y=counts[subset], risktime=risktime[subset], dates=dates[subset]), offset=offset, epsilon=epsilon, max.iter=max.iter)

	

	out        <- list(em=em.model, nlm=nlm.model, misc=misc, measures=results)
	class(out) <- "EMnlm"
	return(out)
}


