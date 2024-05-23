load("compare_methods.results.rda")
ci <- FALSE
#pdf('typical.pdf',width=7,height=7)
par(las=1,bg='white')
if(ci) {
  library(plotrix) # for plotCI
  sx <- sqrt(diag(ret2$posterior$P))
  plotCI(x=truth_condition,y=1+ret2$posterior$x_hat,uiw=sx,err='y',sfrac=0,
         xlab='truth',ylab='retrieved',pch=c(rep('x',264),rep('o',264)),xlim=c(0,2),ylim=c(0,2))
  sx <- sqrt(diag(cov(enkf$dx)))
  plotCI(x=truth_condition,y=enkf$x,uiw=sx,err='y',sfrac=0,add=TRUE,
         col='red',pch=c(rep('x',264),rep('o',264)))
  abline(0,1,col='blue')
  legend(x='topleft',bty='n',pch=c('x','o','x','o',NA),col=c('black','black','red','red','blue'),lty=c(rep(NA,4),1),legend=c('Schuh land','Schuh ocean','EnKF land','EnKF ocean','1:1'))

} else {
#  lims <- range(c(1+ret2$posterior$x_hat,enkf$x,kf$x))
  lims <- range(c(enkf$x,kf$x))
  plot(x=truth_condition,y=1+ret2$posterior$x_hat,xlim=lims,ylim=lims)
  points(x=truth_condition,y=enkf$x,col='red',pch=c(rep('x',264),rep('o',264)))
  points(x=truth_condition,y=kf$x,col='blue',pch=c(rep('x',264),rep('o',264)))
}


#dev.off()
