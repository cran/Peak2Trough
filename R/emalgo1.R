
EMalgo.kfas <- function(
                        ss,
                        Wstruc,			
                        epsilon   = 1e-6,
                        maxiter   = 100,
                        print.ite = FALSE,
                        trace     = FALSE,
                        Wpath     = NULL,
			T         = 2,
			S         = 2,
                        ...
                        ){
	ite <- 0;loglike <- numeric();conv <- 1000
	m0 <- ss$m0;C0 <- ss$C0
	while( conv > epsilon & ite <= maxiter ) {
		ite <- ite + 1
		if(print.ite==TRUE) print(ite)

		# Using KFAS to perform Kalman filtering and smoothing
		smoothed <- sspir::Fkfs(ss, ...)
		R <- smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi) %*% C0 %*% t(smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi)) + smoothed$ss$Wmat(1, smoothed$ss$x, smoothed$ss$phi)
		B <- C0 %*% t(smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi)) %*% solve(R)
		smoothed$ss$C0 <- C0 + B %*% (smoothed$ss$C[[1]] - R) %*% t(B)
		smoothed$ss$m0 <- t(t(m0) + B %*% (smoothed$ss$m[1,] - smoothed$ss$Gmat(1, smoothed$ss$x, smoothed$ss$phi) %*% t(m0)))

		step <- EMstep.kfas(smoothed$ss, m0, C0, T, S, W.est=Wstruc)
		if(ite>1){
			conv <- (smoothed$ss$loglik - loglike[ite-1])
	        }
	
		W <- step$W
		ss$Wmat <- function(tt,x,phi) W


	        if(trace){
			write.table(matrix(step$W[lower.tri(step$W,diag=TRUE)], nrow=1), append = TRUE, col.names = FALSE, row.names = FALSE, file = Wpath)
	        }
		loglike[ite] <- smoothed$ss$loglik
	}

	West <- W
	smoothed$ss$m0 <- m0;  smoothed$ss$C0 <- C0
	fit <- smoothed$ss
	if(conv<epsilon){cc <- TRUE}else{cc <- FALSE}
	list(ss=smoothed$ss, Wmat.est=West, loglik=loglike, iterations=ite, fit=fit, convergence=cc)
}


EMstep.kfas <- function(s, m0, C0, T, S, W.est){
	W.hat <- W.est(s, m0, C0, T, S)
	list(W=W.hat)
}




