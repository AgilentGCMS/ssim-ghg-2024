ct.setup()

# Interesting experiments to try
#
# 1. We don't trust the process model: increase Q to 1
#
# 2. We trust the measurements less than we should. Increase Pd over R.
#
# 3. We trust the measurements too much: Decrease Pd under R.

# scalar, fill rate 0
F <- 1 # state transition, persistence
H <- 1 # measurement mapping: direct to state
L <- runif(min=0,max=100,n=1) # true level
Q <- 0.01 # process noise variance (small)

x0 <- 0 # initial state
P0 <- 1000 # initial state (variance)
R <- 0.06 #  measurement noise
#R <- 0.000006 #  measurement noise
r <- sqrt(R) # std dev of meas noise
n <- 20

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
state$rmse.prior <- rep(NA,n)
state$rmse.post <- rep(NA,n)

nodd <- 0
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
  if(i>1) {
    ylim <- range(c(ylim,c(x[i]+sqrt(Px[i]),x[i]-sqrt(Px[i]))),na.rm=TRUE)
  }
  
  # update

  # noiseless obs
  #d <- L
  d <- L + rnorm(n=1,mean=0,sd=r)
  Pd <- R
  y[i] <- d
  ylim <- range(c(ylim,d),na.rm=TRUE)
  
  # make everything a matrix
  this.x <- matrix(x[i],1,1)
  this.Px <- matrix(Px[i],1,1)
  this.H <- matrix(H,1,1)
  this.d <- matrix(d,1,1)
  this.Pd <- matrix(Pd,1,1)
  xtrue <- matrix(L,1,1)
  post <- measurement.update(this.x,this.Px,this.H,this.d,this.Pd,xtrue)
  x[i] <- post$x[1,1]
  Px[i] <- post$Px[1,1]

  ylim <- range(c(ylim,c(x[i]+sqrt(Px[i]),x[i]-sqrt(Px[i]))),na.rm=TRUE)

  
  state$x.post[i] <- x[i]
  state$Px.post[i] <- Px[i]
  
  state$chi2mn_fcast[i] <- post$chi2mn_fcast
  state$chi2mn_state[i] <- post$chi2mn_state

  state$rmse.post[i] <- compute.rmse(state$x.post[i]-xtrue)
  state$rmse.prior[i] <- compute.rmse(state$x.prior[i]-xtrue)

  if(state$rmse.post[i] > state$rmse.prior[i]) {
    cat(sprintf("[%d] Warning: RMSE post (%.2g) > RMSE prior (%.2g)\n",i,state$rmse.post[i],state$rmse.prior[i]))
    nodd <- nodd+1
  }
  
}

layout(matrix(c(1,1,2,3,4,4),nrow=3,byrow=TRUE))

plot(1:n,y,main="Welch water level model 0: no fill",ylim=ylim,pch='+',xpd=NA)
abline(h=L)
plotCI(x=(1:n)+0.1,y=state$x.post,uiw=sqrt(state$Px.post),col='blue',err='y',add=TRUE,sfrac=0.001)
plotCI(x=(1:n)-0.1,y=state$x.prior,uiw=sqrt(state$Px.prior),col='red',err='y',add=TRUE)

hist(state$chi2mn_fcast,main="Forecast chi square")
abline(v=mean(state$chi2mn_fcast))

hist(state$chi2mn_state,main="State chi square")
abline(v=mean(state$chi2mn_state))

plot(1:n,state$chi2mn_fcast,ylim=range(c(state$chi2mn_fcast,state$chi2mn_state)))
points(1:n,state$chi2mn_state,col='red')
legend(x='topright',bty='n',pch=1,col=c('black','red'),legend=c('fcast','state'))
abline(h=1)

cat(sprintf("chi2mn fcast: %.2g\n",mean(state$chi2mn_fcast)))
cat(sprintf("chi2mn state: %.2g\n",mean(state$chi2mn_state)))
cat(sprintf("no. updates where RMSE is not improved: %d/%d (%.1f%%)\n",
            nodd,n,100*nodd/n))
