ct.setup()

# Water level has a linear trend.

w1 <- function(Q,make.plots=TRUE)  {
  # scalar, fill
  F <- 1 # state transition, persistence
  H <- 1 # measurement mapping: direct to state
  
  #Q <- 0.0001 # process noise variance (small)
  #Q <- 0.0008 # process noise variance (empirically good)
  
  x0 <- 0 # initial state
  P0 <- 1000 # initial state (variance)
  r <- 0.1 # measurement noise (std dev)
  R <- r^2 #  measurement noise (variance)
  n <- 100
  
  y <- rep(NA,n)
  x <- rep(NA,n)
  Px <- rep(NA,n)
  state <- list()
  state$x.prior <- rep(NA,n)
  state$Px.prior <- rep(NA,n)
  state$x.post <- rep(NA,n)
  state$Px.post <- rep(NA,n)
  state$chi2mn_fcast <- rep(NA,n)
  state$chi2mn_state <- rep(NA,n)
  
  # true fill rate
  f <- 0.01
  L <- 1 + f*(1:n)

  ylim <- c(NA,NA)

  for (i in 1:n) {

    # predict
    if(i==1) {
      x[i] <- x0
      Px[i] <- P0
    } else {
      x[i] <- F*x[i-1]
      Px[i] <- F*Px[i-1]*F + Q
    }

    state$x.prior[i] <- x[i]
    state$Px.prior[i] <- Px[i]
    #  if(i>1) {
    #    ylim <- range(c(ylim,c(x[i]+sqrt(Px[i]),x[i]-sqrt(Px[i]))),na.rm=TRUE)
    #  }
    
    # update
    
    d <- L[i] + rnorm(n=1,mean=0,sd=r)
#    d <- L[i]
    Pd <- R
    y[i] <- d
    ylim <- range(c(ylim,d),na.rm=TRUE)
    
    # make everything a matrix
    this.x <- matrix(x[i],1,1)
    this.Px <- matrix(Px[i],1,1)
    this.H <- matrix(H,1,1)
    this.d <- matrix(d,1,1)
    this.Pd <- matrix(Pd,1,1)
    xtrue <- matrix(L[i],1,1)
    post <- measurement.update(this.x,this.Px,this.H,this.d,this.Pd,xtrue)
    x[i] <- post$x[1,1]
    Px[i] <- post$Px[1,1]

    ylim <- range(c(ylim,c(x[i]+sqrt(Px[i]),x[i]-sqrt(Px[i]))),na.rm=TRUE)
    
    state$x.post[i] <- x[i]
    state$Px.post[i] <- Px[i]
    state$chi2mn_fcast[i] <- post$chi2mn_fcast
    state$chi2mn_state[i] <- post$chi2mn_state
    
  }
  
  if(make.plots) {
    layout(matrix(c(1,1,2,3,4,4),nrow=3,byrow=TRUE))

    plot(1:n,y,main="Welch water level model 1: constant fill rate",ylim=ylim)
    points(1:n,L,pch='+')
    plotCI(x=1:n,y=state$x.post,uiw=sqrt(state$Px.post),col='blue',err='y',add=TRUE,sfrac=0.001)
    #plotCI(x=1:n,y=state$x.prior,uiw=sqrt(state$Px.prior),col='red',err='y',add=TRUE)



    hist(state$chi2mn_fcast,main="Forecast chi square")
    abline(v=mean(state$chi2mn_fcast))

    hist(state$chi2mn_state,main="State chi square")
    abline(v=mean(state$chi2mn_state))

    plot(1:n,state$chi2mn_fcast,ylim=range(c(state$chi2mn_fcast,state$chi2mn_state)),log='y')
    points(1:n,state$chi2mn_state,col='red')
    abline(h=1)


    cat(sprintf("Q %.5f - chi2mn fcast: %.2f\n",Q,mean(state$chi2mn_fcast[10:n])))
    cat(sprintf("Q %.5f - chi2mn state: %.2f\n",Q,mean(state$chi2mn_state[10:n])))
  } else {
    retval <- list()
    retval$chi2mn_fcast <- mean(state$chi2mn_fcast[10:n])
    retval$chi2mn_state <- mean(state$chi2mn_state[10:n])
    return(retval)
  }
}

if(TRUE) {
  Q <- 0.0001
  w1(Q=Q)

  invisible(readline(prompt="Press [enter] to continue"))

  Q <- 0.0008
  w1(Q=Q)

} else {

  niter <- 200
  chi2mn_fcast <- list()
  chi2mn_state <- list()

  Qs <- sort(c(0.0007,0.0018,seq(0.0002,0.0010,by=0.0002)))

  x11(width=8,height=20)
  par(las=1,mai=c(0.15,0.1,0.1,0))
  layout(matrix(1:(2*length(Qs)),byrow=TRUE,ncol=2))

  iq <- 0
  for (Q in Qs) {
    iq <- iq+1

    chi2mn_fcast[[iq]] <- rep(NA,niter)
    chi2mn_state[[iq]] <- rep(NA,niter)

    pb <- progress.bar.start(sprintf("Q=%.5f  %d",Q,niter),niter)
    
    for (iter in 1:niter) {
      foo <- w1(Q=Q,make.plots=FALSE)
      chi2mn_fcast[[iq]][iter] <- foo$chi2mn_fcast
      chi2mn_state[[iq]][iter] <- foo$chi2mn_state
      pb <- progress.bar.print(pb,iter)
    }
    progress.bar.end(pb)
    
    hist(chi2mn_fcast[[iq]],main=sprintf("fcast Q=%.5f",Q),breaks=20,col='wheat',xlim=c(0,5))
    abline(v=1,col='black')
    abline(v=mean(chi2mn_fcast[[iq]]),col='brown',lwd=3)
    hist(chi2mn_state[[iq]],main=sprintf("state Q=%.5f",Q),breaks=20,col='wheat',xlim=c(0,5))
    abline(v=1,col='black')
    abline(v=mean(chi2mn_state[[iq]]),col='brown',lwd=3)
  }


}
