## trend - de-season
showtrend <- function(m, origin = 1977, period = 52, conf=TRUE, xlab="", ylab="Incidence rates", titl=NULL, fun = function(m){exp(m)}, ccol = "#3CB37150", ylim=NULL, points=NULL, ...) {

    Axel  <- origin:(m$n/period + origin)
    trend <- m$m[1:m$n,1]
    M <- unlist(lapply( m$C, function(M) M[1,1] ))[1:m$n]
    upper <- trend + 1.96*sqrt(M)
    lower <- trend - 1.96*sqrt(M)

    if(is.null(ylim)){ yy <- range(fun(upper), fun(lower) )}else{ yy <- ylim}
    plot(x=1:length(trend), type="n", main=titl, axes=FALSE, xlab=xlab, ylab=ylab, ylim=yy)
    if(!is.null(points)){
        points(points, pch=20 )
    }

    if(conf){
        polygon(x = c(1:m$n,m$n:1), y = c(fun(lower),rev(fun(upper))), col = ccol, border = NA)
    }

    lines(fun(trend), ...)

    axis(1, at=seq(1,m$n+period,period)[1:length(Axel)], Axel )
    axis(2)
    box()


}


showseason <- function(m, time=1, idx=4:11, origin=1977, titl=NULL, conf=TRUE, period=52, xlab="", ylab="Percentual deviation from annual median", fun = function(m){(exp(m)-1)*100},  ccol = "#3CB37150", ylim = NULL, arrows=FALSE, lcol=1, llty=1, llwd=1) {

    mx <- matrix(m$m[time,idx],ncol=length(idx))
    FF <- matrix(NA,ncol=length(idx),nrow=period)

    for(i in 1:period){

        FF[i,] <- m$Fmat(i, m$x, m$phi)[idx]
    }

    M <- lapply( m$C, function(M) M[idx,idx] )[time]

    res <- matrix(NA, ncol=length(time), nrow=period)
    for(j in 1:length(time)){
        for (i in 1:period){
            res[i,j] <- t(FF[i,]) %*% M[[j]] %*% FF[i,]
        }
    }

    seas  <-  FF %*% t(mx)
    upper <- seas +  1.96*sqrt(res)
    lower <- seas -  1.96*sqrt(res)

    if(is.null(ylim)){ yy <- range(fun(upper), fun(lower) )}else{ yy <- ylim}
    plot(x=1:period, type="n", ylim=yy, main=titl, xlab=xlab, axes=FALSE, ylab=ylab)
    if(conf){
        polygon(c(1:period,period:1), c(fun(lower),rev(fun(upper))), col = ccol, border = NA)
    }
    for(i in 1:length(time)){
        lines(fun(seas[,i]), col=lcol[i], lty=llty[i], lwd=llwd)
    }

    if(arrows){
        maxy <- max(fun(seas))
        maxx <- (1:52)[fun(seas)==maxy]
        miny <- min(fun(seas))
        minx <- (1:52)[fun(seas)==miny]

        arrows(x0=c(maxx, minx), y0=c(min(yy)+3,min(yy)+3), x1=c(maxx, minx), y1=c(min(yy),min(yy)), code=2, col=c("red","blue"), lwd=3)
        arrows(x0=0, y0=maxy, x1=0, y1=miny, code=3, lwd=3, col="green")
    }


    axis(1, at=seq(1,period-1,length.out=12), labels=c(month.abb))
    axis(2)
    box()
}



showlinpred <- function(m, titl=NULL, conf=TRUE, period=52, origin=1977, xlab="", ylab="Incidence rates", fun=function(m){exp(m)},  ccol = "#3CB37150", ylim=NULL, points=NULL, ...){

    Axel  <- origin:(m$n/period + origin)
    FF <- matrix(NA,ncol=m$p,nrow=m$n)
    res <- c()

    for(i in 1:m$n){

        FF[i,] <- m$Fmat(i, m$x, m$phi)
         res <- c(res, t(m$Fmat(i, m$x, m$phi)) %*% m$C[[i]] %*% m$Fmat(i, m$x, m$phi))
    }

    seas  <- diag( FF %*% t(m$m) )
    upper <- seas +  1.96*sqrt(res)
    lower <- seas -  1.96*sqrt(res)

    if(is.null(ylim)){ yy <- range(fun(upper), fun(lower) )}else{ yy <- ylim}

    plot(x=1:length(seas), type="n", ylim=yy, main=titl, xlab=xlab, axes=FALSE, ylab=ylab)

    if(!is.null(points)){
        points(points, pch=20 )
    }

    if(conf){
        polygon(c(1:m$n,m$n:1), c(fun(lower),rev(fun(upper))), col = ccol, border = NA)
    }

    lines(fun(seas), ...)

    axis(1, at=seq(0,m$n+period,period)[1:length(Axel)], Axel )
    axis(2)
    box()

}

"seasonVideo.EMnlm" <- function(ss, ylim, ...){
	animation::saveVideo(
		for(i in 1:ss$nlm$ss$n){
			showseason(ss$nlm$ss, time=i, idx=ss$misc$T+1:(2*ss$misc$S), arrows=TRUE, conf=FALSE, ylim=ylim, titl=ss$misc$data$dates[i])
		}, ...)
}



"plot.EMnlm" <- function(x, ...){

	par(ask=TRUE)
	
	# Plotting 6 random seasonal variation curves
	times <- sort(sample(1:x$nlm$ss$n, 6))
	# Plotting 100*(1 - seasonal variation) (percentual deviation from annual median)
	showseason(x$nlm$ss, time=times, idx=(x$misc$T+1):dim(x$nlm$ss$Wmat(1, x$nlm$ss$x, x$nlm$ss$phi))[2], conf=FALSE, ccol="blue", titl="Seasonal variation", lcol=1:6, llty=1:6, llwd=2, arrows=FALSE)
	legend(x=0, y=0, xjust=0, yjust=1, legend=x$misc$data$dates[times], lty=1:6, col=1:6, lwd=2)

	origin <- as.POSIXlt(as.Date(min(as.numeric(as.Date(x$misc$data$dates, "%d %b %Y"))), origin="1970-01-01"))$year+1900
	
	# Plotting exp(trend)
	showtrend(x$nlm$ss, conf=FALSE, fun=function(m){ exp(m) }, lwd=2,
ylab="Incidence rates", titl="Trend", col="red", ccol="blue", origin=origin, ylim=c(0,max(x$nlm$ss$y/x$misc$data$risktime*x$misc$offset)), points = x$nlm$ss$y/x$misc$data$risktime*x$misc$offset)

	# Plotting exp(linear predictor)
	showlinpred(x$nlm$ss, conf=FALSE, fun=function(m){ exp(m) }, lwd=2, origin=origin, ylab="Incidence rates", titl="Linear predictor", ylim=c(0,max(x$nlm$ss$y/x$misc$data$risktime*x$misc$offset)), points = x$nlm$ss$y/x$misc$data$risktime*x$misc$offset, col="red", ccol="blue")
}

