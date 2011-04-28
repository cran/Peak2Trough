rr1 <-
function (counts, plot = TRUE) 
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
    if (plot) {
        plot(counts, xlab = "Time", ylab = "Count")
        li <- sum(counts)/12 * (1 + alpha * cos(((0.5:((5 * k) - 
            0.5)) * 2 * pi/(5 * k)) - thetamax))
        lines(y = li, x = (0.5:((5 * k) - 0.5))/5)
        points(y = max(li), x = Month, pch = 19)
    }
    list(RelativeRisk = rr, TimePeak = Month)
}

