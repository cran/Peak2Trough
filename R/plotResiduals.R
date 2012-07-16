"plot.residuals" <- function(
	resid, 
	path  = NULL, 
	mtextx=c(expression(T[t-1]), expression(a[t-1]),
		expression(alpha[list(1,t-1)]), expression(beta[list(1,t-1)]),
		expression(alpha[list(2,t-1)]), expression(beta[list(2,t-1)]),
		expression(tilde(T)[t-1]), expression(tilde(a)[t-1]),
		expression(tilde(alpha)[list(1,t-1)]), expression(tilde(beta)[list(1,t-1)]),
		expression(tilde(alpha)[list(2,t-1)]), expression(tilde(beta)[list(2,t-1)])), 
	mtexty=c(expression(T[t]), expression(a[t]), expression(alpha[list(1,t)]),
		expression(beta[list(1,t)]), expression(alpha[list(2,t)]),
		expression(beta[list(2,t)]), expression(tilde(T)[t]),
		expression(tilde(a)[t]), expression(tilde(alpha)[list(1,t)]),
		expression(tilde(beta)[list(1,t)]), expression(tilde(alpha)[list(2,t)]),
		expression(tilde(beta)[list(2,t)])), 
	burn.in=1, ...){

	n         <- resid$misc$ss$n
	p         <- resid$misc$ss$p
	FyRes     <- resid$FyRes
	SyRes     <- resid$SyRes
	FthetaRes <- resid$FthetaRes
	SthetaRes <- resid$SthetaRes
	FITy      <- resid$misc$FITy 
	FITtheta  <- resid$misc$FITtheta
	FITsy     <- resid$misc$FITsy 
	FITstheta <- resid$misc$FITstheta
	

        # Markov assumption origin y
        pdf(file=paste(path, '/v lag1.pdf', sep=""), ...)
        plot(y=FyRes[(burn.in+1):n], x=FyRes[burn.in:(n-1)], xlab="", ylab="", pch=20)
        mtext(side=1, line=2, expression(v[t-1]))
        mtext(side=2, expression(v[t]), line=2)
        dev.off()
        pdf(file=paste(path, '/vtilde lag1.pdf', sep=""), ...)
        plot(y=SyRes[(burn.in+1):n], x=SyRes[burn.in:(n-1)], xlab="", ylab="", pch=20)
        mtext(side=1, line=2, expression(tilde(v)[t-1]))
        mtext(side=2, expression(tilde(v)[t]), line=2)
        dev.off()

        # Markov assumption origin theta
        for(i in 1:p){
		pdf(file=paste(path, '/w', i, 'lag1.pdf', sep=""), ...)
		plot(y=FthetaRes[i,(burn.in+1):n], x=FthetaRes[i,burn.in:(n-1)], xlab="", ylab="", pch=20)
		mtext(side=1, line=2, mtextx[i])
		mtext(side=2, mtexty[i], line=2)
		dev.off()
		pdf(file=paste(path, '/wtilde', i, 'lag1.pdf', sep=""), ...)
		plot(y=SthetaRes[i,(burn.in+1):n], x=SthetaRes[i,burn.in:(n-1)], xlab="", ylab="", pch=20)
		mtext(side=1, line=2, mtextx[i+p])
		mtext(side=2, mtexty[i+p], line=2)
		dev.off()
	}

        # Distribution assumption origin y
        pdf(file=paste(path, '/v dist.pdf', sep=""), ...)
        plot(y=FyRes[burn.in:n], x=FITy[burn.in:n], xlab="", ylab="", pch=20)
        mtext(side=1, line=2, expression(paste(F[t]^T,G[t],m[t-1])))
        mtext(side=2, line=2, expression(v[t]))
        dev.off()
        pdf(file=paste(path, '/vtilde dist.pdf', sep=""), ...)
        plot(y=SyRes[burn.in:n], x=FITsy[burn.in:n], xlab="", ylab="", pch=20)
        mtext(side=1, line=2, expression(paste(F[t]^T,tilde(m)[t])))
        mtext(side=2, line=2, expression(tilde(v)[t]))
        dev.off()

        # Distribution assumption origin theta
        for(i in 1:p){
		pdf(file=paste(path, '/w', i, 'dist.pdf', sep=""), ...)
		plot(y=FthetaRes[i,burn.in:n], x=FITtheta[i,burn.in:n], xlab="", ylab="", pch=20)
		mtext(side=1, line=2, expression(paste(G[t],m[t-1])))
		mtext(side=2, line=2, mtexty[i])
		dev.off()
		pdf(file=paste(path, '/wtilde', i, 'dist.pdf', sep=""), ...)
		plot(y=SthetaRes[i,burn.in:n], x=FITstheta[i,burn.in:n], xlab="", ylab="", pch=20)
		mtext(side=1, line=2, expression(paste(G[t],tilde(m)[t-1])))
		mtext(side=2, line=2, mtexty[i+p])
		dev.off()
        }

        # ACF plots
        pdf(file=paste(path, '/v acf.pdf', sep=""), ...)
        acf(FyRes[burn.in:n], lag.max=52, xlab="", main="")
        mtext(side=1, line=2, expression(v[t]))
        dev.off()
        pdf(file=paste(path, '/vtilde acf.pdf', sep=""), ...)
        acf(SyRes[burn.in:n], lag.max=52, xlab="", main="")
        mtext(side=1, line=2, expression(tilde(v)[t]))
        dev.off()

        for(i in 1:p){
		pdf(file=paste(path, '/w', i, 'acf.pdf', sep=""), ...)
		acf(FthetaRes[i,burn.in:n], lag.max=52, xlab="", main="")
		mtext(side=1, line=2, mtexty[i])
		dev.off()
		pdf(file=paste(path, '/wtilde', i, 'acf.pdf', sep=""), ...)
		acf(SthetaRes[i,burn.in:n], lag.max=52, xlab="", main="")
		mtext(side=1, line=2, mtexty[i+p])
		dev.off()
        }

        # Time plots
        pdf(file=paste(path, '/v time.pdf', sep=""), ...)
        plot(y=FyRes[burn.in:n], x=burn.in:n, xlab="Time", ylab="", pch=20, cex=1)
        mtext(side=2, line=3, expression(v[t]))
        dev.off()
        pdf(file=paste(path, '/vtilde time.pdf', sep=""), ...)
        plot(y=SyRes[burn.in:n], x=burn.in:n, xlab="Time", ylab="", pch=20, cex=1)
        mtext(side=2, line=3, expression(tilde(v)[t]))
        dev.off()

        for(i in 1:p){
		pdf(file=paste(path, '/w', i, 'time.pdf', sep=""), ...)
		plot(y=FthetaRes[i,burn.in:n], x=burn.in:n, xlab="Time", ylab="", pch=20, cex=1)
		mtext(side=2, line=2, mtexty[i])
		dev.off()
		pdf(file=paste(path, '/wtilde', i, 'time.pdf', sep=""), ...)
		plot(y=SthetaRes[i,burn.in:n], x=burn.in:n, xlab="Time", ylab="", pch=20, cex=1)
		mtext(side=2, line=2, mtexty[i+p])
		dev.off()
        }
}

