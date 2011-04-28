#Simulating data
RR    <- c(1,1.10,1.5,2)
sample.size <- c(150,500,2500)
num.ite <- 10
# In order to reconstruct results from paper set num.ite <- 1000
                     
dta <- dta.pool <- list()
Bias.ed.1 <- mse.ed.1 <- Bias.poi.1 <- mse.poi.1 <- Bias.poifull.1 <- mse.poifull.1 <- 
matrix(c(0,150,500,2500,1.00,rep(0,3),1.10,rep(0,3),1.50,rep(0,3),2.00,rep(0,3)),nrow=4)

set.seed(678)
k <- 1
for (i in 1:3){
for (j in 1:4){
temp <- model1(RR=RR[j],years=28,num.dta=num.ite,sample.size=sample.size[i])
dta <- temp$sim.data
dta.pool <- temp$sim.agr.data
RR.ed <- RR.poi <- RR.poi.f <- matrix(NA,ncol=num.ite,nrow=24)
  for(l in 1:num.ite){
    RR.ed[k,l] <- rr1(dta.pool[[l]],plot=FALSE)$RelativeRisk
    RR.poi[k,l] <- rr2(dta[[l]],plot=FALSE)$RelativeRisk
    RR.poi.f[k,l] <- rr3(dta[[l]],plot=FALSE)$RelativeRisk
  }
  Bias.ed.1[(i+1),(j+1)] <- mean(RR.ed[k,]) - Bias.ed.1[1,(j+1)]
  mse.ed.1[(i+1),(j+1)]  <- sd(RR.ed[k,])
  Bias.poi.1[(i+1),(j+1)] <- mean(RR.poi[k,]) - Bias.poi.1[1,(j+1)]
  mse.poi.1[(i+1),(j+1)]  <- sd( RR.poi[k,] )
  Bias.poifull.1[(i+1),(j+1)] <- mean(RR.poi.f[k,]) - Bias.poifull.1[1,(j+1)]
  mse.poifull.1[(i+1),(j+1)]  <- sd( RR.poi.f[k,] )
k <- k+1
}                                             
}

N150   <- formatC(cbind(Bias.ed.1[1,],100*Bias.ed.1[2,],100*mse.ed.1[2,],100*Bias.poi.1[2,],100*mse.poi.1[2,],100*Bias.poifull.1[2,],100*mse.poifull.1[2,])[2:5,], format="f", digits = 2, drop0trailing=FALSE,flag=" ")
N500   <- formatC(cbind(Bias.ed.1[1,],100*Bias.ed.1[3,],100*mse.ed.1[3,],100*Bias.poi.1[3,],100*mse.poi.1[3,],100*Bias.poifull.1[3,],100*mse.poifull.1[3,])[2:5,], format="f", digits = 2, drop0trailing=FALSE,flag=" ")
N2500  <- formatC(cbind(Bias.ed.1[1,],100*Bias.ed.1[4,],100*mse.ed.1[4,],100*Bias.poi.1[4,],100*mse.poi.1[4,],100*Bias.poifull.1[4,],100*mse.poifull.1[4,])[2:5,], format="f", digits = 2, drop0trailing=FALSE,flag=" ")

#N=150
N150
#N=500
N500
#N=2500
N2500
