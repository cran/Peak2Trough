"residuals.EMnlm" <- function(ss){

	m0  <- ss$nlm$ss$m0
	C0  <- ss$nlm$ss$C0
	fm  <- ss$nlm$ss$mf
	fC  <- ss$nlm$ss$Cf
	sm  <- ss$nlm$ss$m
	sC  <- ss$nlm$ss$C
	G   <- ss$nlm$ss$Gmat
	F   <- ss$nlm$ss$Fmat
	W   <- ss$nlm$ss$Wmat(1, ss$nlm$ss$x, ss$nlm$ss$phi)
	x   <- ss$nlm$ss$x
	phi <- ss$nlm$ss$phi
	V   <- ss$nlm$ss$vtilde
	y   <- ss$nlm$ss$ytilde

	R   <- as.matrix(nearPD(G(1,x,phi) %*% C0 %*% t(G(1,x,phi)) + W)$mat) #t
	B   <- C0 %*% t(G(1,x,phi)) %*% solve( R )             #t-1
	sC0 <- as.matrix(nearPD(C0 + B %*% (sC[[1]] - R) %*% t(B))$mat)
	sm0 <- t(t(m0) + B %*% (sm[1,] - G(1, x, phi) %*% t(m0)))
	L   <- as.matrix(nearPD(sC[[1]] + G(1,x,phi) %*% sC0 %*% t(G(1,x,phi)) - sC[[1]] %*% t( B ) %*% t(G(1,x,phi)) - G(1,x,phi) %*% B %*% t(sC[[1]]))$mat)

	Q   <- t(F(1,x,phi)) %*% R %*% F(1,x,phi) + V[1]
U <- as.matrix(nearPD(R %*% matrix(F(1,x,phi), ncol=1) %*% solve(Q) %*% t(matrix(F(1,x,phi), ncol=1)) %*% R)$mat)

	v   <- sqrt(Q)^(-1) * (y[1] - t(F(1,x,phi)) %*% G(1,x,phi) %*% matrix(m0))
	w   <- msqi(U) %*% (fm[1,] - G(1,x,phi) %*% matrix(m0))
	Q   <- V[1] - t(F(1,x,phi)) %*% sC[[1]] %*% F(1,x,phi)
	U   <- as.matrix(nearPD(W - L)$mat)

	v.tilde <- sqrt(Q)^(-1) * (y[1] - t(F(1,x,phi)) %*% sm[1,])
	w.tilde <- msqi(U) %*% (sm[1,] - G(1,x,phi) %*% matrix(sm0))
	
	fit.y       <- t(F(1,x,phi)) %*% G(1,x,phi) %*% matrix(m0)
        fit.theta   <-   G(1,x,phi)  %*% matrix(m0)
        fit.smo     <- t(F(1,x,phi)) %*% sm[1,]
        fit.s.theta <-   G(1,x,phi)  %*% matrix(sm0)

	for(t in 2:ss$nlm$ss$n){
	        R <- as.matrix(nearPD(G(t,x,phi) %*% fC[[t-1]] %*% t(G(t,x,phi)) + W)$mat) #t
	        B <- fC[[t-1]] %*% t(G(t,x,phi)) %*% solve(R)  #t-1
	        L <- as.matrix(nearPD(sC[[t]]+ G(t,x,phi) %*% sC[[t-1]] %*% t(G(t,x,phi)) - sC[[t]] %*% t(B) %*% t(G(t,x,phi)) - G(t,x,phi) %*% B %*% t(sC[[t]]))$mat)

	        Q <- t(F(t,x,phi)) %*% R %*% F(t,x,phi) + V[t]
	        U <- as.matrix(nearPD(R %*% F(t,x,phi) %*% solve(Q) %*% t(F(t,x,phi)) %*% R)$mat)

	        v <- c(v, sqrt(Q)^(-1) * (y[t] - t(F(t,x,phi)) %*% G(t,x,phi) %*% fm[t-1,]))
	        w <- cbind(w, msqi(U) %*% (fm[t,] - G(t,x,phi) %*% fm[t-1,]))
	        Q <- V[t] - t(F(t,x,phi)) %*% sC[[t]] %*% F(t,x,phi)
	        U <- as.matrix(nearPD(W - L)$mat)

	        v.tilde <- c(v.tilde, sqrt(Q)^(-1) * (y[t] - t(F(t,x,phi)) %*% sm[t,]))
	        w.tilde <- cbind(w.tilde, msqi(U) %*% (sm[t,] - G(t,x,phi) %*% sm[t-1,]))

		fit.y <- c(fit.y, t(F(t,x,phi)) %*% G(1,x,phi) %*% fm[t-1,])
		fit.theta <- cbind(fit.theta, G(t,x,phi) %*% fm[t-1,])
		fit.smo <- c(fit.smo, t(F(t,x,phi)) %*% sm[t,])
		fit.s.theta <- cbind(fit.s.theta, G(t,x,phi) %*% sm[t-1,])
	}

	misc=list(ss=ss$nlm$ss, FITy=fit.y, FITtheta=fit.theta, FITsy=fit.smo, FITstheta=fit.s.theta)
	
	res=list(FyRes=v, FthetaRes=w, SyRes=v.tilde, SthetaRes=w.tilde, misc=misc)
	class(res) <- "residuals"
	return(res)
}

#msqi <- function(m){
#
#    eig <- eigen(m)
#    isq <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% solve(eig$vectors)
#    return(isq)
#
#}

msqi <- function(m){
    diag(1/sqrt(diag(m)))
}

#msqi <- function(m){
# diag(10)
#}
