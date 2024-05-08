ctn0sfr.time.stamp <- "Time-stamp: <tropics.local:/Users/andy/Presentations/2024/enkf_summer_school/notebooks/welch/kf_support_functions.r - 07 Mar 2024 (Thu) 09:32:34 MST>"

compute.rmse <- function(x,na.rm=TRUE) {
  return(sqrt(mean(x^2,na.rm=na.rm)))
}


# Kalman filter measurement update.
#
# update the solution {x, Sx} with additional observations {z, Sz}
# where measurements are related to state via H %*% x = z
#
# If xtrue is supplied, compute chi2mn_state
#
# Allows for measurement rejection if threshold is supplied. That
# threshold is in number of standard deviations.

measurement.update <- function(x,Sx,H,z,Sz,xtrue=NULL,
                               rejection.threshold=NA,
                               verbose=FALSE) {
  
  err <- FALSE
  msgs <- character(0)
  
  if(is.null(dim(z))) {
    msgs <- c(msgs,sprintf("observations array z has no dimensions, needs to be nobs x 1 matrix"))
    err <- TRUE
  }
  if(is.null(dim(H))) {
    msgs <- c(msgs,sprintf("H matrix has no dimensions, needs to be nobs x nparms matrix"))
    err <- TRUE
  }

  nparms <- dim(x)[1]
  nobs <- dim(z)[1]

  if(!all(dim(Sx)==nparms)) {
    msgs <- c(msgs,sprintf("nparms is %d, but dim(Sx) is %d,%d",nparms,dim(Sx)[1],dim(Sx)[2]))
    err <- TRUE
  }

  if(!all(dim(x)==c(nparms,1))) {
    msgs <- c(msgs,sprintf("nparms is %d, but dim(x) is %d,%d",nparms,dim(x)[1],dim(x)[2]))
    err <- TRUE
  }

  if(!all(dim(z)==c(nobs,1))) {
    msgs <- c(msgs,sprintf("nobs is %d, but dim(d) is %d,%d",nobs,dim(z)[1],dim(z)[2]))
    err <- TRUE
  }
  
  if(!all(dim(H)==c(nobs,nparms))) {
    msgs <- c(msgs,sprintf("nparms is %d, nobs is %d but dim(H) is %d,%d",nparms,nobs,dim(H)[1],dim(H)[2]))
    err <- TRUE
  }

  if(!all(dim(Sz)==nobs)) {
    msgs <- c(msgs,sprintf("nobs is %d, but dim(Sz) is %d,%d",nobs,dim(Sz)[1],dim(Sz)[2]))
    err <- TRUE
  }
  
  if(err) {
    for (imsg in 1:length(msgs)) {
      cat(sprintf("[measurement.update] %s\n",msgs[imsg]))
    }
    stop("Error(s) in measurement.update()")
  }

  # if xtrue, save prior x minus xtrue as prior x error
  x.pri <- x
  Sx.pri <- Sx
  if(!is.null(xtrue)) {
    ex.pri <- x-xtrue
  }

  retval <- list()
  retval$assimilated <- rep(TRUE,nobs)

  retval$fcast_error <- (z - H %*% x)
  Fm <-  H %*% Sx %*% t(H) + Sz
  Fmic <- solve(t(chol(Fm)))
  retval$fcast_error_normalized <- Fmic %*% retval$fcast_error
  nobs.orig <- nobs
  
  if(is.na(rejection.threshold)) {
    nobs.rej <- 0
    jobs.assim <- nobs
    lx.ass <- 1:nobs
    lx.rej <- numeric(0)
  } else {
    Fm <-  H %*% Sx %*% t(H) + Sz
    Fmi <- solve(Fm)

    lx.rej <- which(abs(retval$fcast_error_normalized)>rejection.threshold)
    lx.ass <- which(abs(retval$fcast_error_normalized)<=rejection.threshold)
    retval$assimilated[lx.rej] <- FALSE
    nobs.rej <- length(lx.rej)
    nobs.assim <- length(lx.ass)
    z <- z[lx.ass]
    Sz <- Sz[lx.ass,lx.ass]
    H <- H[lx.ass,]
    nobs <- length(lx.ass)
    cat(sprintf("[measurement.update] Rejecting %d of %d obs (%.1f%%)\n",
                nobs.rej,nobs.orig,100*nobs.rej/nobs.orig))
  }

  # Now post-rejection quantities
  
  fcast_error <- (z - H %*% x)
  Fm <-  H %*% Sx %*% t(H) + Sz
  Fmi <- solve(Fm)
  K <- (Sx %*% t(H)) %*% Fmi
  deltax <- K %*% fcast_error
  x <- x + deltax
  retval$x <- x
  imkh <- diag(rep(1,length(x)))-K %*% H

  # Both of these posterior Sx forms are supposed to be equivalent
  
  #retval$Sx <- imkh %*% Sx %*% t(imkh) + K %*% Sz %*% t(K)
  Sx <- imkh %*% Sx
  retval$Sx <- Sx

  retval$kappa <- kappa(retval$Sx)

  # could add diagnostics right here for normality of parameters
  # mean(x) should be zero
  # distribution of x should be multivariate normal N(0,Sx)
  
  # Compute no. of DOFs following Patil et al (2001)
  if(TRUE){
    ndofs <- nparms-1
  } else {
    library(svd)
    dd <- svd(Sx)
    ndofs <- round(sum(dd$d)^2/sum(dd$d^2))
  }
  retval$ndofs <- ndofs
  retval$chi2mn_fcast <- (1/nobs) * t(fcast_error) %*% Fmi %*% fcast_error
  retval$post_error <- rep(NA,nobs.orig)
  retval$post_error[lx.ass] <- (z - H %*% x )

  if(verbose) {
      if(compute.rmse(retval$post_error) > compute.rmse(retval$fcast_error)) {
          stop("Measurement update worsened forecast errors\n")
      }
  }

  if(sum(diag(Sx)) > sum(diag(Sx.pri))) {
    stop("Trace of posterior covariance exceeds prior")
  }
  
  if(!is.null(xtrue)) {
      ex <- x-xtrue
      retval$x.error <- ex
    Sxi <- solve(Sx) #inverse of sqrt of Sx prior for chi2 analysis only
    retval$chi2mn_state <- (1/ndofs) * t(ex) %*% Sxi %*% ex

    if(verbose) {
        if((t(ex) %*% ex) > (t(ex.pri) %*% ex.pri)) {
            cat("Measurement update worsened agreement with xtrue\n")
        }
    }

  }
  
  return(retval)
  
}

#  Ensenble Kalman filter measurement update.
#
# update the solution x (with uncertainty ensemble dx)
# with additional observations d
enkf.measurement.update <- function(x,dx,d,Pd,y,dy,xtrue=NULL) {

  # In Whitaker & Hamill 2002 (henceforth, whitaker02a),
  # their H x'b is our dy,
  # their H Pb Ht is our HQH (Q is covariance matrix before assimilation),
  # their x'b is our dx  (parameter deviations before assimilation).

  # NOTE R's chol() returns an upper-triangular matrix, whereas
  # python's cholesky() returns a lower-triangular matrix.  Looks
  # like lower-triangular is the correct result, so all R calls to
  # chol() are wrapped in t().

  err <- FALSE
  msgs <- character(0)

  nmemb <- dim(dy)[2]
  nparms <- dim(x)[1]
  nobs <- dim(d)[1]

  if(!all(dim(x)==c(nparms,1))) {
    msgs <- c(msgs,sprintf("nparms is %d, but dim(x) is %d,%d",nparms,dim(x)[1],dim(x)[2]))
    err <- TRUE
  }
  
  # check dimensions
  
  #  dy nobs, nmemb
  if(!all(dim(dy)==c(nobs, nmemb))) {
    msgs <- c(msgs,sprintf("dim(dy) expected to be nobs, nmemb (%d, %d), but got (%d, %d)",
                           nobs,nmemb,dim(dy)[1],dim(dy)[2]))
    err <- TRUE
  }
  
  #  x nparms, 1
  if(!all(dim(x)==c(nparms, 1))) {
    msgs <- c(msgs,sprintf("dim(x) expected to be x nparms, 1 (%d, %d), but got (%d, %d)",
                               nparms,1,dim(x)[1],dim(x)[2]))
              err <- TRUE
    }
    
  #  dx nparms, nmemb
  if(!all(dim(dx)==c(nparms, nmemb))) {
    msgs <- c(msgs,sprintf("dim(dx) expected to be nparms, nmemb (%d, %d), but got (%d, %d)",
                           nparms,nmemb,dim(dx)[1],dim(dx)[2]))
    err <- TRUE
  }

  #  Pd nobs, nobs
  if(!all(dim(Pd)==c(nobs, nobs))) {
    msgs <- c(msgs,sprintf("dim(Pd) expected to be  nobs, nobs (%d, %d), but got (%d, %d)",
                               nobs,nobs,dim(Pd)[1],dim(Pd)[2]))
    err <- TRUE
  }

  #  d nobs,1
  if(!all(dim(d)==c(nobs,1))) {
    msgs <- c(msgs,sprintf("dim(d) expected to be nobs,1 (%d, %d), but got (%d, %d)",
                                nobs,1,dim(d)[1],dim(d)[2]))
    err <- TRUE
  }
  
  HQHR <- (1/(nmemb-1))*(dy %*% t(dy)) + Pd
    
  if(err) {
    for (imsg in 1:length(msgs)) {
      cat(sprintf("[measurement.update] %s\n",msgs[imsg]))
    }
    stop("Error(s) in enkf.measurement.update()")
  }

  # if xtrue, save prior x minus xtrue as prior x error
  x.pri <- x
  #  Px.pri <- Px
  if(!is.null(xtrue)) {
    ex.pri <- x-xtrue
  }

  HQHRi <- solve(HQHR)
  HQ <- (1/(nmemb-1))*(dx %*% t(dy))

  # invert central value
  
  K <- HQ %*% HQHRi
  
  # apply localization mask to K
#  K <- localization_mask * K # element-by-element multiplication
 
  resid <- d-y  # obs minus sim, y being the simulated obs from the "mean" member
  x <- x + K %*% resid

  sHQHR <- t(chol(HQHR))
#  sHQHRi <- ginv(sHQHR)
  sHQHRi <- solve(sHQHR)

  Ktilde <- HQ %*% t(sHQHRi) %*% solve(sHQHR + t(chol(Pd)))  # whitaker02a, eqn 10

  # apply localization mask to Ktilde
  #Ktilde <- localization_mask * Ktilde # element-by-element multiplication
  
  dx <- dx - Ktilde %*% dy

  retval <- list()
  retval$dx <- dx

  retval$x <- x
  
  return(retval)
  
}
