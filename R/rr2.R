rr2 <-
function (counts, period = 365, offz = NULL, plot = TRUE, origin="1960-01-01", hemisphere=c("Northern","Southern"), x.axis=c("Months","Seasons"), y.range=list(NULL,NULL))
{
    k <- period
    angles <- (1:length(counts)) * 2 * pi/k
    if(length(offz)==0) { offz <- rep(1,length(counts)) }
    temp.model <- glm(counts ~ cos(angles) + sin(angles), x = TRUE,
        family = poisson(), offset=log(offz))
    rr <- exp(2 * sqrt(temp.model$coefficient[2]^2 + temp.model$coefficient[3]^2))
    rd <- exp(sqrt(temp.model$coefficient[2]^2 + temp.model$coefficient[3]^2)) -
          exp(-sqrt(temp.model$coefficient[2]^2 + temp.model$coefficient[3]^2))
    thetamax <- atan(abs(temp.model$coefficient[3]/temp.model$coefficient[2]))
    thetamax <- sign(temp.model$coefficient[3]) * sign(temp.model$coefficient[2]) *
        thetamax
    thetamax <- pi * (temp.model$coefficient[2] < 0) + (temp.model$coefficient[2] >
        0) * (temp.model$coefficient[3] < 0) * 2 * pi + thetamax
    month <- (thetamax * k)/(2 * pi)
    if (plot) {
      if (hemisphere=="Northern") {
        ## Northern hemisphere ###
              par(mfrow=c(2,1))
              plot(counts/offz, xlab = "Time", ylab = "Incidence rate",pch = "*", main = "Exponentiated linear predictor", ylim=c(y.range[[1]][1],y.range[[1]][2]) , axes=FALSE)
              axis(1,seq(1,length(counts)+k,k), c(format(as.Date(seq(1,length(counts)+k,k)*365/k, origin=origin),"%Y")))
              axis(2)
              box()
              li <- exp(temp.model$x %*% temp.model$coefficient)
              lines(y = li , x = (1:length(counts)), lwd = 3, col = "red")
          if(x.axis=="Months"){
              plot(y = li[1:k], x = 1:k,  xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[2]][1],y.range[[2]][2]), axes=FALSE)
              axis(1,cumsum(rep(k/12,13))-k/12,c(month.abb,""))
              axis(2)
              box()
              points(y = max(li), x = month, pch = 19 , col="red")
          } else {
              lii <- c(li[(k-round(k/12)+1):k],li[1:(k-round(k/12))])
              plot(y = lii , x = 1:k,  xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[2]][1],y.range[[2]][2]), axes=FALSE)
              axis(1,seq(0,k-round(k/4),k/4), c("Winter", "Spring", "Summer", "Autumn"),tck=0)
              axis(2)
              box()
              points(y = max(li), x = (1:k)[lii=max(lii)], pch = 19 , col="red")
          }
      } else {

       ## Southern hemisphere ###
        par(mfrow=c(2,1))
        plot(counts/offz, xlab = "Time", ylab = "Incidence rate",pch = "*", main = "Exponentiated linear predictor", ylim=c(y.range[[1]][1],y.range[[1]][2]) , axes=FALSE)
        axis(1,seq(1,length(counts)+k,k),
             c(format(as.Date(seq(1,length(counts)+k,k)*365/k,
             origin=origin),"%Y")))
        axis(2)
        box()
        lines(y = li , x = (1:length(counts)), lwd = 3, col = "red")
        if(x.axis=="Months"){
            lii <-  c(li[(floor(k/12*6)):k],li[1:(floor(k/12*6)-1)])
            plot(lii , x = 1:k,  xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[2]][1],y.range[[2]][2]), axes=FALSE)
            axis(1,cumsum(rep(k/12,13))-k/12, c(month.abb[7:12], month.abb[1:6],""))
            axis(2)
            box()
            points(y = max(li), x = (1:k)[lii==max(lii)], pch = 19 , col="red")
        } else {
            li <-  c(li[(floor(k/12*6)):k],li[1:(floor(k/12*6)-1)])
            lii <- c(li[(k-round(k/12)+1):k],li[1:(k-round(k/12))])
            plot(y = lii , x = 1:k, xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[2]][1],y.range[[2]][2]), axes=FALSE)
            axis(1,seq(0,k-round(k/4),k/4)+k/8, c("Winter", "Spring", "Summer", "Autumn"),tck=0)
            axis(2)
            box()
            points(y = max(li), x = (1:k)[lii==max(lii)], pch = 19 , col="red")
      }
    }
  }

    cat("------------------------------- \n")
    cat("Relative risk   = \t",  formatC(rr,format="f", digits = 2, drop0trailing=FALSE,flag="#",width=6) , "\n" )
    cat("Risk difference = \t",  formatC(rd,format="f", digits = 2, drop0trailing=FALSE,flag="#",width=6) , "\n" )
    cat("Peaktime        =\t", formatC(format(as.Date(365/k*month, origin="1960-01-01"),"%d/%m"),format="s", flag=" ",width=6) , "\n")
    cat("------------------------------- \n")

    list(Model=temp.model, RelativeRisk = rr, RiskDif=rd,  TimePeak = format(as.Date(365/k*month, origin="1960-01-01"),"%d/%m"))
}

