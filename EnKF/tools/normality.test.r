normality.test <- function(e,plot=TRUE,verbose=TRUE,known.mean=NULL,known.sd=NULL) {

  lx <- which(is.na(e))
  if(length(lx)>0) {
    
    n <- length(e)
    if(verbose) {
      cat(sprintf("\nWARNING: Removing %d missing values.\n",length(lx)))
    }
    e <- e[setdiff(1:n,lx)]
  }
  
  n <- length(e)
  if(plot) {
    layout(matrix(1:4,2,2,byrow=T))
    plot(e,type='l',col='grey')
    points(e,pch=21,cex=0.5,bg="black")
    hist(e,min((n/2),50))
    qqnorm(e)
    qqline(e,col='red')
  }
  
  
  # Apply a battery of tests to determine whether or not the
  # given vector is plausibly a set of independent, identically
  # disctributed samples of a normal Gaussian population

  retval <- list()
  retval$e <- e
  
  if(verbose) {
    cat('\nH0: vector e is normal; H1: e is non-normal.\n')
  }
  retval$H0 <- "e is normal"
  retval$H1 <- "e is not normal"
  
  library(stats)
  if(n<5000) {
    retval$shapiro <- shapiro.test(e)
  } else {
    retval$shapiro <- list()
    retval$shapiro$text <- "too many residuals for test"
    retval$shapiro$p.value <- NA
  }
  if(!is.null(known.mean)&!is.null(known.sd)) {
    retval$kolmogorov <- ks.test(e,"pnorm",mean=known.mean,sd=known.sd)
  }

  # REFERENCES FOR SKEWNESS AND KURTOSIS
  #
  #  Box 7.1, section 7.9 and box 7.5 of:
  #
  # @Book{sokal95a,
  #   author =	 {Robert R. Sokal and F. James Rohlf},
  #   title = 	 {Biometry},
  #   publisher = 	 {Freeman},
  #   address =	 {New York},
  #   year = 	 1995,
  #   edition =	 {Third}
  # }
  #
  #   Sokal and Rohlf reference D'Agostino and collaborators:
  #

  # @Article{dagostino73a,
  #   author = 	 {Ralph B. {D'Agostino} and Gary L. Tietjen},
  #   title = 	 {Approaches to the null distribution of $\sqrt{b_1}$},
  #   journal = 	 {Biometrika},
  #   year = 	 1973,
  #   volume =	 60,
  #   pages =	 {169--173}
  # }
  # 
  # @Article{dagostino73b,
  #   author = 	 {Ralph B. {D'Agostino} and Gary L. Tietjen},
  #   title = 	 {Tests for departure from normality. {Empirical} results for the distribution of $b_2$ and $\sqrt{b_1}$},
  #   journal = 	 {Biometrika},
  #   year = 	 1973,
  #   volume =	 60,
  #   pages =	 {613-622}
  # }
  # 

  # Report mean and variance
  if(verbose) {
    cat('N:',n,' Mean:',mean(e),'Variance:',var(e),'\n')
  }
  retval$n <- n
  retval$mean <- mean(e)
  retval$var <- var(e)
  

  # Assess skewness
  library(e1071)
  sk <- skewness(e)
  retval$skew <- sk
  if(verbose) {
    cat('Skewness:',sk)
    if(sk<0) {
      cat(' (left skew)')
    } else if(sk>0) {
      cat(' (right skew)')
    }
  }
  if(sk<0) {
    retval$skew.text <- 'left skew'
  } else if(sk>0) {
    retval$skew.text <- 'right skew'
  }
  
  sg1 <- sqrt((6*n*(n-1))/((n-2)*(n+1)*(n+3)));  # standard error of the skewness (g1)

  t.sk <- (sk - 0)/sg1;

  # significance values for this t statistic come from the t CDF for
  # infinite DOF.  For infinite DOF, t converges to the normal distribution.

  pval <- 2*(1-pnorm(abs(t.sk)));
  if(verbose) {
    cat(" p=",formatC(pval,digits=4,format='f'),"\n")
  }
  retval$p.skew <- pval
  
  # Assess kurtosis

  kt <- kurtosis(e);
  retval$kurtosis <- kt
  
  if(verbose) {
    cat('Kurtosis',kt)
    if(kt<0) {
      cat(' (platykurtic)')
    } else if (kt>0) {
      cat(' (leptokurtic)')
    }
  }
  if(kt<0) {
    retval$kurtosis.text <- "platykurtic"
  } else if (kt>0) {
    retval$kurtosis.text <- "leptokurtic"
  }
  
  sg2 <- sqrt((24*n*(n-1)*(n-1))/((n-3)*(n-2)*(n+3)*(n+5)));  # standard error of the skewness (g2)

  t.kt <- (kt - 0)/sg2;

  # significance values for this t statistic come from the t CDF for
  # infinite DOF.  For infinite DOF, t converges to the normal distribution.

  pval = 2*(1-pnorm(abs(t.kt)));

  if(verbose) {
    cat(" p=",formatC(pval,digits=4,format='f'),"\n")
  }
  retval$p.kurt <- pval
  
  library(stats)
  eacf <- acf(e,lag.max=n-2,plot=plot,main="ACF")
  retval$eacf <- eacf

  # Tong:  no more than 5% of ACF members beyond lag zero should exceed +/- 1.96/sqrt(n)
  acf.thresh <- 1.96/sqrt(n)
  n.exceed <- length(which(abs(eacf$acf[2:n])>acf.thresh))
  pct.exceed <- 100*n.exceed/(n-1)
  retval$acf.n.exceed <- n.exceed
  retval$acf.pct.exceed <- pct.exceed
  
  invisible(retval)
  
}
