ct.setup()

# Water level has a linear trend.
# state is now two-dimensional: level and fill rate

# measurement mapping: state is observed, fill rate unobserved
H <- matrix(c(1,0),nrow=1)

# process noise variance
#Q <- 0.015*matrix(c(0.01,0,0,0.01),ncol=2) 
Q <- 24*matrix(c(0.008,0,0,0.013),ncol=2) 
#Q <- matrix(c(0.0001,0,0,0.000005),ncol=2) 

x0 <- matrix(c(0,0),ncol=1) # initial state
P0 <- matrix(c(1000,0,0,1000),ncol=2) # initial state (variance)

r <- 1 # measurement noise (std dev)
R <- r^2 #  measurement noise (variance)
n <- 20
y <- rep(NA,n) # obs

state <- list()
state$x.prior <- array(NA,dim=c(2,n))
state$Px.prior <- array(NA,dim=c(2,2,n))
state$x.post <- array(NA,dim=c(2,n))
state$Px.post <- array(NA,dim=c(2,2,n))
state$chi2mn_fcast <- rep(NA,n)
state$chi2mn_state <- rep(NA,n)
state$rmse_post <- rep(NA,n)
state$rmse_prior <- rep(NA,n)

# true fill rate
#f <- 0.5
f <- rnorm(mean=0.5,sd=0.1,n=1)
#f <- rnorm(mean=0.5,sd=0.1,n=n)
state$xtrue <- array(NA,dim=c(2,n))
state$xtrue[1,] <- 1 + f*(1:n) # state (level)
state$xtrue[2,] <- rep(f,n) # fill rate
#state$xtrue[1,] <- 1 + f # state (level)
#tate$xtrue[2,] <- f

deltat <- 1
psi <- matrix(c(1,0,deltat,1),ncol=2) # persistence except x1[i+1] <- x1[i] + x2[i]*deltat

nodd <- 0
ylims <- list()
ylims[[1]] <- c(NA,NA)
ylims[[2]] <- c(NA,NA)

for (i in 1:n) {

  # predict
  if(i==1) {
    state$x.prior[,i] <- x0
    state$Px.prior[,,i] <- P0
  } else {
    state$x.prior[,i] <- psi %*% state$x.post[,i-1]
    state$Px.prior[,,i] <- psi %*% state$Px.post[,,i-1] %*% t(psi) + Q
  }

  # update
  
  d <- H %*% state$xtrue[,i] + rnorm(n=1,mean=0,sd=r)
#  d <- H %*% state$xtrue[,i]
  Pd <- R
  y[i] <- d
  ylims[[1]] <- range(c(ylims[[1]],d),na.rm=TRUE)
  
  # make everything a matrix
  this.d <- matrix(d,1,1)
  this.Pd <- matrix(Pd,1,1)

  post <- measurement.update(matrix(state$x.prior[,i],nrow=2,ncol=1),
                             state$Px.prior[,,i],
                             H,
                             this.d,this.Pd,state$xtrue[,i])
  state$x.post[,i] <- post$x
  state$Px.post[,,i] <- post$Px

  state$chi2mn_fcast[i] <- post$chi2mn_fcast
  state$chi2mn_state[i] <- post$chi2mn_state

  state$rmse.post[i] <- compute.rmse(state$x.post[,i]-state$xtrue[,i])
  state$rmse.prior[i] <- compute.rmse(state$x.prior[,i]-state$xtrue[,i])

  if(state$rmse.post[i] > state$rmse.prior[i]) {
    cat(sprintf("[%d] Warning: RMSE post (%.2g) > RMSE prior (%.2g)\n",i,state$rmse.post[i],state$rmse.prior[i]))
    nodd <- nodd+1
  }
  ylims[[1]] <- range(c(ylims[[1]],
                        c(post$x[1]+sqrt(post$Px[1,1]),post$x[1]-sqrt(post$Px[1,1]))),
                      na.rm=TRUE)

  ylims[[2]] <- range(c(ylims[[2]],
                        c(post$x[2]+sqrt(post$Px[2,2]),post$x[2]-sqrt(post$Px[2,2]))),
                      na.rm=TRUE)
  
  
}

par(mai=c(0.5,0.3,0.1,0.1))
layout(matrix(c(1,1,2,2,3,4,5,5),nrow=4,byrow=TRUE))

plot(1:n,y,main="Welch water level model 2: x[1]=level",ylim=ylims[[1]])
points((1:n),state$xtrue[1,],pch='+')
plotCI(x=(1:n)+0.2,y=state$x.post[1,],uiw=sqrt(state$Px.post[1,1,]),col='blue',err='y',add=TRUE,sfrac=0.001)
plotCI(x=(1:n)-0.2,y=state$x.prior[1,],uiw=sqrt(state$Px.prior[1,1,]),col='red',err='y',add=TRUE)
legend(x='topright',bty='n',col=c('black','red','blue',"black"),legend=c('true','prior','post','obs'),pch=c(3,1,1))

plotCI(x=2:n,y=state$x.post[2,2:n],uiw=sqrt(state$Px.post[2,2,2:n]),col='blue',err='y',sfrac=0.001,main="Welch water level model 2: x[2]=fill rate",xlim=c(1,n),ylim=c(-1,2),xpd=NA)
#points(x=3:n,f[3:n],col='black')
abline(h=f)


hist(state$chi2mn_fcast,main="Forecast chi square")
abline(v=c(1,mean(state$chi2mn_fcast)),col=c('black','red'))

hist(state$chi2mn_state,main="State chi square")
abline(v=c(1,mean(state$chi2mn_state)),col=c('black','red'))

plot(1:n,state$chi2mn_fcast,ylim=range(c(state$chi2mn_fcast,state$chi2mn_state)),log='y')
points(1:n,state$chi2mn_state,col='red')
legend(x='topright',bty='n',pch=1,col=c('black','red'),legend=c('fcast','state'))
abline(h=1)

cat(sprintf("chi2mn fcast: %.2f\n",mean(state$chi2mn_fcast)))
cat(sprintf("chi2mn state: %.2f\n",mean(state$chi2mn_state)))
cat(sprintf("no. updates where RMSE is not improved: %d/%d (%.1f%%)\n",
            nodd,n,100*nodd/n))
