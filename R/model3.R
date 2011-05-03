model3 <-
function (RR = 1, years = 1, num.dta = 10, daily = TRUE, sample.size = 500, 
    Xs, Bs) 
{
    if (daily) {
        k <- 365
    }
    else {
        k <- 12
    }
    dato <- (1:(k * years))
    hh <- (Xs %*% Bs)[1:k]
    f <- (log(RR))/(max(hh) - min(hh))  
    prob <- list()
    temp <- list()
    temp.pool <- list()
    for (l in 1:num.dta) {
        prob[[l]] <- exp((Xs %*% Bs) * f)
        temp[[l]] <- rmultinom(1, sample.size * years, prob[[l]])
        dta.temp <- data.frame(count = temp[[l]])
        dta.temp$month <- format(as.Date(dato, origin = "1960-01-01"), 
            "%m")
        temp.pool[[l]] <- as.matrix(aggregate(dta.temp$count, 
            by = list(dta.temp$month), FUN = sum)[, 2], nrow = 12)
    }
    list(sim.data = temp, sim.agr.data = temp.pool)
}

