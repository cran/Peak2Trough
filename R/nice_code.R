W2phi <- function(W, T, S){
	t.phi <- t.W2phi(W[1:T,1:T])
	s.phi <- s.W2phi(W[(T+1):dim(W)[2],(T+1):dim(W)[2]])
	c(t.phi,s.phi)
}

phi2x <- function(phi, T, S){
	W <- phi2W(phi, T, S)
	if(W[dim(W)[2],dim(W)[2]]==0){ 
		cor          <- matrix(0, ncol=dim(W)[2], nrow=dim(W)[1])
		cor[1:T,1:T] <- cov2cor(W[1:T,1:T])
		diag(cor) <- c(log(diag(W[1:T,1:T])), rep(0,2*S))
	}else{
		cor <- cov2cor(W)
		diag(cor) <- log(diag(W))
	}
	cor[lower.tri(cor)] <- logit((cor[lower.tri(cor)]+1)/2)
	cor[upper.tri(cor)] <- logit((cor[upper.tri(cor)]+1)/2)
	W2phi(cor, T, S)
}

x2phi <- function(x, ...){
	cor <- phi2W(x, ...)
	sd <- sqrt(exp(diag(cor)))
	diag(cor) <- 1
	cor[lower.tri(cor)] <- 2*expit(cor[lower.tri(cor)])-1
	cor[upper.tri(cor)] <- 2*expit(cor[upper.tri(cor)])-1
	W <- outer(sd,sd)*cor
	W2phi(W,...)
}

phi2W <- function(phi, T, S){
	t.W <- t.phi2W(phi[1:(T*(T+1)/2)])
	s.W <- s.phi2W(tail(phi, n=2*S+S*(S-1)/2))
	W <- matrix(0, T+2*S,T+2*S)
	W[1:T,1:T] <- t.W
	W[(T+1):dim(W)[2],(T+1):dim(W)[2]] <- s.W
	W
}

H.trans <- function(T, S){
	t.H <- diag(T)
	s.H <- rbind(kronecker(diag(S),matrix(c(1, 1),1,2)),kronecker(diag(S),matrix(c(1, -1),1,2)))
	H <- matrix(0, T+2*S,T+2*S)
	H[1:T,1:T] <- t.H
	H[(T+1):dim(H)[2],(T+1):dim(H)[2]] <- s.H
	H
}




t.phi2W <- function(phi){

t <- (sqrt(1+8*length(phi))-1)/2
W <- matrix(0, t, t)
W[lower.tri(W)] <- phi[(t+1):length(phi)]

W <- W+t(W)

diag(W) <- phi[1:t]

W
}


s.phi2W <- function(phi){

d <- (sqrt(9+8*length(phi))-3)/2
Wsmall <- matrix(0, d, d)
Wsmall[lower.tri(Wsmall)] <- phi[(2*d+1):length(phi)]
Wsmall <- Wsmall+t(Wsmall)
diag(Wsmall) <- phi[(d+1):(2*d)]
W <- kronecker(Wsmall,matrix(1,2,2))
diag(W) <- kronecker(phi[1:d],matrix(1,2,1))

W

}


t.W2phi <- function(W){

t <- dim(W)[1]
vectorized <- as.numeric(W)

index <- numeric(t*(t-1)/2)
for(i in 1:t){
	temp <- 1:i
	index[i*(i-1)/2+temp] <- temp+(i-1)*t
}

c(diag(W), vectorized[-index])

}

s.W2phi <- function(W){

d <- dim(W)[1]/2
Wreduce <- W[-seq(1,2*d,by=2),-seq(2,2*d,by=2)]
temp <- t.W2phi(Wreduce)
Wdiag <- diag(W)

c(Wdiag[seq(1,2*d,by=2)],temp)

}


W.strc <- function(s, m0, C0, T, S){
    x <- s$x
    phi <- s$phi

    R <- s$Gmat(1,x,phi) %*% C0 %*% t(s$Gmat(1,x,phi)) + s$Wmat(1,x,phi) #t
    B <- C0 %*% t(s$Gmat(1,x,phi)) %*% solve( R )             #t-1
    L <- s$C[[1]]+ s$Gmat(1,x,phi) %*% s$C0 %*% t(s$Gmat(1,x,phi)) -
       s$C[[1]] %*% t( B ) %*% t(s$Gmat(1,x,phi)) -
       s$Gmat(1,x,phi) %*% B %*% t(s$C[[1]])

    W1 <- L+(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1))) %*%
        t(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1)))


    for(t in 2:(s$n)){
        R <- s$Gmat(t,x,phi) %*% s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) + s$Wmat(1,x,phi) #t
        B <- s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) %*% solve( R )  #t-1
        L <- s$C[[t]]+ s$Gmat(t,x,phi) %*% s$C[[t-1]] %*% t(s$Gmat(t,x,phi)) - s$C[[t]] %*% t( B ) %*% t(s$Gmat(t,x,phi)) - s$Gmat(t,x,phi) %*% B %*% t(s$C[[t]])

        W1 <- W1 + L + (t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1))) %*% t(t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1)))

    }

    W <- matrix(0, nrow=T+2*S, ncol=T+2*S)
    W[1:T,1:T] <- W1[1:T,1:T]
    W[(T+1):(T+S),(T+1):(T+S)] <- W1[(T+1):(T+S),(T+1):(T+S)]
    diag(W[(T+1+S):dim(W)[2],(T+1+S):dim(W)[2]]) <- diag(W1[(T+1+S):dim(W)[2],(T+1+S):dim(W)[2]])

    W <- 1/(2*s$n)*(W + t(W))

    return(W)

}

W.strc.stat <- function(s, m0, C0, T, S){
	x   <- s$x
	phi <- s$phi

	R   <- s$Gmat(1,x,phi) %*% C0 %*% t(s$Gmat(1,x,phi)) + s$Wmat(1,x,phi) #t
	B   <- C0 %*% t(s$Gmat(1,x,phi)) %*% solve( R )             #t-1
	L   <- s$C[[1]]+ s$Gmat(1,x,phi) %*% s$C0 %*% t(s$Gmat(1,x,phi)) - s$C[[1]] %*% t( B ) %*% t(s$Gmat(1,x,phi)) - s$Gmat(1,x,phi) %*% B %*% t(s$C[[1]])

	W1  <- L+(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1))) %*% t(t(matrix(s$m[1,],nrow=1)) - s$Gmat(1,x,phi) %*% t(matrix(s$m0[1,],nrow=1)))


	for(t in 2:(s$n)){
        	R  <- s$Gmat(t,x,phi) %*% s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) + s$Wmat(1,x,phi) #t
        	B  <- s$Cf[[t-1]] %*% t(s$Gmat(t,x,phi)) %*% solve( R )  #t-1
        	L  <- s$C[[t]]+ s$Gmat(t,x,phi) %*% s$C[[t-1]] %*% t(s$Gmat(t,x,phi)) - s$C[[t]] %*% t( B ) %*% t(s$Gmat(t,x,phi)) - s$Gmat(t,x,phi) %*% B %*% t(s$C[[t]])
        	W1 <- W1 +L+(t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1))) %*% t(t(matrix(s$m[t,], nrow=1)) - s$Gmat(t,x,phi) %*% t(matrix(s$m[t-1,], nrow=1)))
	}

	W   <- matrix(0, nrow=T+2*S, ncol=T+2*S)
	W[1:T,1:T] <- W1[1:T,1:T]
    
	W   <- 1/(2*s$n)*(W + t(W))

	return(W)
}




smooth.bottom <- function(ss){
	R <- as.matrix(nearPD(ss$Gmat(1, ss$x, ss$phi) %*% ss$C0 %*% t(ss$Gmat(1, ss$x, ss$phi)) + ss$Wmat(1, ss$x, ss$phi))$mat) #t
	B <- ss$C0 %*% t(ss$Gmat(1, ss$x, ss$phi)) %*% solve( R )             #t-1
	sC0 <- as.matrix(nearPD(ss$C0 + B %*% (ss$C[[1]] - R) %*% t(B))$mat)
	sm0 <- t(t(ss$m0) + B %*% (ss$m[1,] - ss$Gmat(1, ss$x, ss$phi) %*% t(ss$m0)))

	list(m0=sm0, C0=sC0)
}

sspir2kfas <- function(M, dim, ...){
    ## Tiny help function
    Mt <- vapply(1:dim[3], M, array(0,dim=c(dim[1],dim[2])), ...)
}

calc.season <- function(ss, time=1, idx=3:6, period=52, fun = function(m){exp(m)}) {

	mx <- matrix(ss$m[time,idx],ncol=length(idx))
	xx <- seq(1,52,length.out=52*7)
	FF <- matrix(NA,ncol=4,nrow=length(xx))

	for(i in 1:length(xx)){
	        FF[i,] <- c(cos(2*pi*xx[i]/period) , sin(2*pi*xx[i]/period) , cos(4*pi*xx[i]/period) , sin(4*pi*xx[i]/period) )
	    }


	seas     <-  fun(FF %*% t(mx))
	time.min <- c(1:length(xx))[seas==min(seas)]
	time.max <- c(1:length(xx))[seas==max(seas)]
	ptt      <- max(seas)/min(seas)
	ptt.dif  <- max(seas)-min(seas)

	names(ptt) <- "PTR";names(ptt.dif) <- "PTD"
	c(ptt, ptt.dif, time.min, time.max)
}
