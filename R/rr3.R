rr3 <-
function (counts, daily = TRUE, plot = TRUE) 
{
    if (daily) {
        k <- 365
    }
    else {
        k <- 12
    }
    angles <- (1:length(counts)) * 2 * pi/k
    temp.model <- glm(counts ~ rcs(c(1:length(counts)), 5) + 
        cos(angles) + sin(angles) + cos(angles * 2) + sin(angles * 
        2) + cos(angles * 3) + sin(angles * 3) + cos(angles * 
        4) + sin(angles * 4), x = TRUE, family = poisson())
    g <- exp(temp.model$x[, c(6:13)] %*% matrix(temp.model$coefficients)[c(6:13)])[1:k]
    rr <- max(g)/min(g)
    month <- c(1:k)[g == max(g)]
    if (plot) {
        par(mfrow = c(3, 1))
        plot(counts, xlab = "Time", ylab = "Count", pch = "*", 
            main = "Exponentiated linear predictor")
        li <- exp(temp.model$x %*% temp.model$coefficient)
        lines(y = li, x = (1:length(counts)), lwd = 3, col = "red")
        plot(exp(temp.model$x[,1:5] %*% matrix(temp.model$coefficient[1:5])),
            xlab = "Time", ylab = "Count", lwd = 3, 
            main = "Exponentiated secular trend")
        plot(y = (g - 1) * 100, x = 1:365, type = "l", xlab = "Time", 
            ylab = "Percentual deviation from median", lwd = 3, 
            main = "Detrended seasonal variation")
        points(y = max((g - 1) * 100), x = month, pch = 19, col = "red")
    }
    list(RelativeRisk = rr, TimePeak = month)
}
