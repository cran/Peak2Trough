rr3 <-
function (counts, dates, period = 365, offz = NULL, plot = TRUE, origin="1960-01-01", hemisphere=c("Northern","Southern"), x.axis=c("Months","Seasons"), y.range=list(NULL,NULL,NULL))
{
    k <- period
    if(length(offz)==0) { offz <- rep(1,length(counts)) }

    temp.model <- glm(counts ~ rcs(dates, 5) +
                               cos(2*pi*dates/k) + sin(2*pi*dates/k) + 
                               cos(2*pi*dates/k * 2) + sin(2*pi*dates/k * 2) + 
                               cos(2*pi*dates/k * 3) + sin(2*pi*dates/k * 3) +  
                               cos(2*pi*dates/k * 4) + sin(2*pi*dates/k * 4), 
                      x = TRUE, family = poisson(), offset=offz)

    g <- exp(temp.model$x[, c(6:13)] %*% matrix(temp.model$coefficients)[c(6:13)])[1:k]
    rr <- max(g)/min(g)
    rd <- max(g)-min(g)
    month.max <- c(1:k)[g == max(g)]
    month.min <- c(1:k)[g == min(g)]
    if (plot) {
      if (hemisphere=="Northern") {
          ## Northern hemisphere ###
          par(mfrow = c(3, 1))
          plot(counts/offz, xlab = "Time", ylab = "Incidence rate", pch = "*", main = "Exponentiated linear predictor", ylim = c(y.range[[1]][1], y.range[[1]][2]), axes=FALSE)
          axis(1,seq(1,length(counts)+k,k), c(format(as.Date(seq(1,length(counts)+k,k)*365/k, origin=origin),"%Y")))
          axis(2)
          box()
          li <- exp(temp.model$x %*% temp.model$coefficient)
          lines(y = li, x = (1:length(counts)), lwd = 3, col = "red")
          plot(exp(temp.model$x[,1:5] %*% matrix(temp.model$coefficient[1:5])), xlab = "Time", ylab = "Incidence rate", lwd = 3, type="l", main = "Exponentiated secular trend", ylim = c(y.range[[2]][1], y.range[[2]][2]), axes=FALSE)
          axis(1,seq(1,length(counts)+k,k),
             c(format(as.Date(seq(1,length(counts)+k,k)*365/k,
             origin=origin),"%Y")))
          axis(2)
          box()
          li <- (g-1)*100
          if(x.axis=="Months"){
              plot(y = li , x = 1:k, type = "l", xlab = "", ylab = "Percentual deviation from median", lwd = 3, ylim = c(y.range[[3]][1], y.range[[3]][2]), axes=FALSE, main = "Detrended seasonal variation")
              axis(1,cumsum(rep(k/12,13))-k/12,c(month.abb,""))
              axis(2)
              box()
              points(y = max((g - 1) * 100), x = month.max, pch = 19, col = "red")
          } else {
              lii <- c(li[(k-round(k/12)+1):k],li[1:(k-round(k/12))])
              plot(y = lii, x = 1:k,  xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[3]][1],y.range[[3]][2]), axes=FALSE)
              axis(1,seq(0,k-round(k/4),k/4)+k/8, c("Winter", "Spring", "Summer", "Autumn"),tck=0)
              axis(2)
              box()
              points(y = max(li), x = (1:k)[lii==max(lii)], pch = 19 , col="red")
          }

      } else {

        ## Southern hemisphere ###
        par(mfrow = c(3, 1))
        plot(counts/offz, xlab = "Time", ylab = "Incidence rate", pch = "*", main = "Exponentiated linear predictor", ylim = c(y.range[[1]][1], y.range[[1]][2]), axes=FALSE)
        axis(1,seq(1,length(counts)+k,k),
             c(format(as.Date(seq(1,length(counts)+k,k)*365/k,
             origin=origin),"%Y")))
        axis(2)
        box()
        li <- exp(temp.model$x %*% temp.model$coefficient)
        lines(y = li, x = (1:length(counts)), lwd = 3, col = "red")
        plot(exp(temp.model$x[,1:5] %*% matrix(temp.model$coefficient[1:5])),
            xlab = "Time", ylab = "Incidence rate", lwd = 3, type="l",
            main = "Exponentiated secular trend", ylim = c(y.range[[2]][1], y.range[[2]][2]), axes=FALSE)
        axis(1,seq(1,length(counts)+k,k),
             c(format(as.Date(seq(1,length(counts)+k,k)*365/k,
             origin=origin),"%Y")))
        axis(2)
        box()
        li <- (g - 1) * 100
        if(x.axis=="Months"){
            plot(li , x = 1:k,  xlab = "", ylab = "Percentual deviation from median", lwd=3, type="l", main = "Detrended seasonal variation", ylim=c(y.range[[3]][1],y.range[[3]][2]), axes=FALSE)
            axis(1,cumsum(rep(k/12,13))-k/12, c(month.abb[7:12], month.abb[1:6],""))
            axis(2)
            box()
            points(y = max(li), x = (month.max-(k/12*6))%%k, pch = 19 , col="red")
        } else {
            li <-  c(li[(floor(k/12*6)):k],li[1:(floor(k/12*6)-1)])
            lii <- c(li[(k-round(k/12)+1):k],li[1:(k-round(k/12))])
            plot(y = lii , x = 1:k, xlab = "", ylab = "Count", lwd=3, type="l", main="Exponentiated seasonal variation", ylim=c(y.range[[3]][1],y.range[[3]][2]), axes=FALSE)
            axis(1,seq(0,k-round(k/4),k/4)+k/8, c("Winter", "Spring", "Summer", "Autumn"),tck=0)
            axis(2)
            box()
            points(y = max(li), x =(1:k)[lii==max(lii)] , pch = 19 , col="red")
        }
    }
  }

    end <- list(Model=temp.model, RelativeRisk = rr, RiskDif=rd, TimePeak = format(as.Date(365/k*month.max, origin="1960-01-01"),"%d/%m"), Troughtime = formatC(format(as.Date(365/k*month.min, origin="1960-01-01"),"%d/%m"),format="s", flag=" ",width=6))
    class(end) <- "rr"
    end
}
