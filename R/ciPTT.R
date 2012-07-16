bootstrap.res <- function(fit, nboot=1000){
	b.res <- matrix(sample(fit$Model$residuals, nboot, replace=TRUE), ncol=nboot, nrow=length(fit$Model$residuals))
	fit$Model$y + b.res
}

bootstrap.block <- function(y, b.length=3, nboot=1000){
	options(warn=-1)
	n <- length(y)
	l <- n-b.length+1
	k <- round(n/b.length)
	v <- c(rep(TRUE, b.length), rep(FALSE, l))
	M <- matrix(v, nrow=l, ncol=n, byrow=TRUE)
	blocks <- y * t(M)
	blocks[blocks==matrix(0, ncol=l, nrow=n)] <- NA
	boot <- matrix(sample(1:l, k*nboot, replace=TRUE), ncol=nboot, nrow=k)
	options(warn=0)
	matrix(c(blocks[,c(boot)])[!is.na(c(blocks[,c(boot)]))], nrow=k*b.length, ncol=nboot)
}


"ciPTT" <- function(y, nboot, b.length=10, estimator, alpha=c(0.025, 0.975)){
	n <- length(y)
	thetahat <- estimator(y)
	new.y <- bootstrap.block(y, b.length, nboot)

	est <- numeric()
	for(i in 1:nboot){
		est[i] <- estimator(new.y[,i])
	}

	z0 <- qnorm(sum(est < thetahat)/nboot)
	u <- rep(0, n)
	for (i in 1:n) {
		u[i] <- estimator(y[-i])
	}
	uu <- mean(u) - u
	acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
	zalpha <- qnorm(alpha)
	tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
	ooo <- trunc(tt * nboot)
	confpoints <- sort(est)[ooo]
	confpoints <- cbind(alpha, confpoints)
	dimnames(confpoints)[[2]] <- c("alpha", "bca point")
	return(list(confpoints = confpoints, z0 = z0, acc = acc, u = u))

}

"ciPTT.EMnlm" <- function(m, nboot, b.length=10, estimator, alpha=c(0.025, 0.975)){
	thetahat <- numeric()

	for(i in 1:m$nlm$ss$n){
		thetahat[i] <- calc.season(m$nlm$ss, time=i, period=m$misc$period, idx=m$misc$T+(1:(m$misc$S*2)))[1]
	}
	
	new.y <- bootstrap.block(c(m$nlm$ss$y), b.length, nboot)

	est <- matrix(NA, ncol=m$nlm$ss$m, nrow=nboot)
	for(i in 1:nboot){
		m$nlm$ss$y <- matrix(new.y[,i], ncol=1)
		
		fit <- sspir::Fkfs(m$nlm$ss, tvar=c(m$nlm$ss$n, 1,1,1), offset=m$misc$offset)
		est[i,] <- calc.season(m$nlm$ss, time=i, period=m$misc$period, idx=m$misc$T+(1:(m$misc$S*2)))[1]
	}

	z0 <- qnorm(sum(est < thetahat)/nboot)
	u <- rep(0, m$nlm$ss$n)
	for (i in 1:m$nlm$ss$n) {
		u[i] <- estimator(m$nlm$ss$y[-i])
	}
	uu <- mean(u) - u
	acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
	zalpha <- qnorm(alpha)
	tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
	ooo <- trunc(tt * nboot)
	confpoints <- sort(est)[ooo]
	confpoints <- cbind(alpha, confpoints)
	dimnames(confpoints)[[2]] <- c("alpha", "bca point")
	return(list(confpoints = confpoints, z0 = z0, acc = acc, u = u))

}


