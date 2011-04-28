rr2 <-
function (counts, daily = TRUE, plot = TRUE) 
{
    if (daily) {
        k <- 365
    }
    else {
        k <- 12
    }
    angles <- (1:length(counts)) * 2 * pi/k
    temp.model <- glm(counts ~ cos(angles) + sin(angles), x = TRUE, 
        family = poisson())
    rr <- exp(2 * sqrt(temp.model$coefficient[2]^2 + temp.model$coefficient[3]^2))
    thetamax <- atan(abs(temp.model$coefficient[3]/temp.model$coefficient[2]))
    thetamax <- sign(temp.model$coefficient[3]) * sign(temp.model$coefficient[2]) * 
        thetamax
    thetamax <- pi * (temp.model$coefficient[2] < 0) + (temp.model$coefficient[2] > 
        0) * (temp.model$coefficient[3] < 0) * 2 * pi + thetamax
    month <- (thetamax * k)/(2 * pi)
    if (plot) {
        plot(counts, xlab = "Time", ylab = "Count")
        li <- exp(temp.model$x %*% temp.model$coefficient)
        lines(y = li, x = 1:k)
        points(y = max(li), x = month, pch = 19)
    }
    list(RelativeRisk = rr, TimePeak = month)
}

