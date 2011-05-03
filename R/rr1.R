rr1 <-
function (counts, plot = TRUE, hemisphere=c("Northern","Southern"), x.axis=c("Months","Seasons"), y.range=NULL)
{
    k <- length(as.vector(counts))
    if (k != 12)
        stop("Data is not (aggregated) monhtly counts")
    N <- sum(counts)
    midtheta <- (0.5:(k - 0.5)) * 2 * pi/k
    C <- sum(cos(midtheta) * counts)/k
    S <- sum(sin(midtheta) * counts)/k
    D <- sqrt(C^2 + S^2)
    f <- ((D^2 * k^2)/N)/(1 + (D^2 * k^2)/N)
    thetamax <- atan(abs(S/C))
    thetamax <- sign(S) * sign(C) * thetamax
    thetamax <- pi * (C < 0) + (C > 0) * (S < 0) * 2 * pi + thetamax
    Month <- (thetamax * k)/(2 * pi)
    alpha <- 2 * sqrt((D^2 * k^2 - N * f)/(N * (N - 1)))
    rr <- (1 + alpha)/(1 - alpha)
    rd <- (1 + alpha)-(1 - alpha)
    if (plot) {
      if (hemisphere=="Northern") {
          if(x.axis=="Months"){
             ## Northern hemisphere ###
             plot(y=counts, x=0:11, xlab = "", ylab = "Count", ylim=c(y.range[1],y.range[2]), axes=FALSE)
             axis(1,0:11,month.abb)
             axis(2)
             box()
             li <- sum(counts)/12 * (1 + alpha * cos(seq(0, 11, length.out=60) * 2 * pi/12 - thetamax))
             lines(y = li, x = seq(0, 11, length.out=60), lwd = 3)
             points(y = max(li), x = Month, pch = 19 , col="red")
          } else {
             plot(y=c(counts[12],counts[1:11]), x=0:11, xlab = "", ylab = "Count", ylim=c(y.range[1],y.range[2]), axes=FALSE)
             axis(1, c(1,4,7,10), c("Winter","Spring","Summer","Autumn"),tck=0)
             axis(2)
             box()
             li <- sum(counts)/12 * (1 + alpha * cos(seq(0, 11, length.out=60) * 2 * pi/12 - (thetamax+2*pi/12)))
             lines(y = li, x = seq(0, 11, length.out=60), lwd = 3)
             points(y = max(li), x = Month+1, pch = 19 , col="red")
          }

      } else {
        ## Southern hemisphere ###
          if(x.axis=="Months"){
              plot(y=c(counts[7:12],counts[1:6]), x=0:11, xlab = "", ylab = "Count", ylim=c(y.range[1],y.range[2]), axes=FALSE)
              axis(1,0:11,c(month.abb[7:12],month.abb[1:6]))
              axis(2)
              box()
              li <- sum(counts)/12 * (1 + alpha * cos(seq(0, 11, length.out=60) * 2 * pi/12 - (thetamax+pi)))
              lines(y = li, x = seq(0, 11, length.out=60), lwd = 3)
              points(y = max(li), x = Month%%6, pch = 19 , col="red")
          } else {
             plot(y=c(counts[6:12],counts[1:5]), x=0:11, xlab = "", ylab = "Count", ylim=c(y.range[1],y.range[2]), axes=FALSE)
             axis(1, c(1,4,7,10), c("Winter","Spring","Summer","Autumn"),tck=0)
             axis(2)
             box()
             li <- sum(counts)/12 * (1 + alpha * cos(seq(0, 11, length.out=60) * 2 * pi/12 - (thetamax+7*pi/6)))
             lines(y = li, x = seq(0, 11, length.out=60), lwd = 3)
             points(y = max(li), x = Month%%6+1, pch = 19 , col="red")
         }
      }
  }


    cat("------------------------------- \n")
    cat("Relative risk   = \t",  formatC(rr,format="f", digits = 2, drop0trailing=FALSE,flag="#",width=6) , "\n" )
    cat("Risk difference = \t",  formatC(rd,format="f", digits = 2, drop0trailing=FALSE,flag="#",width=6) , "\n" )
    cat("Peaktime        =\t", formatC(format(as.Date(365/12*Month, origin="1960-01-01"),"%d/%m"),format="s", flag=" ",width=6) , "\n")
    cat("------------------------------- \n")

    list(Data=counts, RelativeRisk = rr, RiskDif=rd, TimePeak = format(as.Date(365/12*Month, origin="1960-01-01"),"%d/%m"))
}

