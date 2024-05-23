# Time-stamp: <hercules-login-1.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/tools/enkf.r: 22 May 2024 (Wed) 22:31:52 UTC>

compute.rmse <- function(x,na.rm=TRUE) {
  return(sqrt(mean(x^2,na.rm=na.rm)))
}

# This function computes a fast Pearson's r
my.Pearsons.corr <- function(dx,dy) {
  n <- dim(dx)[2]

  t1 <- n*tcrossprod(dx,dy)
  t2 <- rowSums(dx) %o% rowSums(dy)
  
  t3 <- sqrt(n*rowSums(dx*dx)-(rowSums(dx))^2)
  t4 <- sqrt(n*rowSums(dy*dy)-(rowSums(dy))^2)
  
  my.rho <- (t1-t2)/(t3 %o% t4) # verified as good

  return(my.rho)
}

ndofs.patil <- function(S) {
# Compute no. of DOFs following Patil et al (2001)
  library(svd)
  dd <- svd(S)
  ndofs <- round(sum(dd$d)^2/sum(dd$d^2))
  return(ndofs)
}

svd.pseudoinverse <- function(x,tol=1e-8) {
  library(svd)
  xd <- svd(x)
  middle.term <- diag(1/xd$d)
  lx <- which(xd$d/xd$d[1] < tol)
  if(length(lx)>0) {
    cat(sprintf("[svd.pseudoinverse] Truncating %d SVs\n",length(lx)))
    middle.term[lx,lx] <- 0
    stop()
  }
  return(xd$v%*%middle.term%*%t(xd$u))
}

generate_ensemble <- function(x=NA,Sx,nmemb) {

  # Compute ensemble of random samples from the multivariate
  # distribution specified by {x,Sx}. Returns a matrix of
  # nmemb,nparms. Note this can also be done via MASS' mvrnorm().
  
  # x is the (nparms vector) mean of the desired distribution
  #   (if x not provided, assume zero)
  # Sx is the (nparms x nparms) covariance matrix,
  # nmemb is the number of samples

  # R's chol() returns an upper-triangular matrix; many applications
  # assume that the result is lower-triangular. t() not applied as it
  # doesn't matter for this application.
  D <- chol(Sx)
  
  nparms <- dim(Sx)[1]

  # generate sequence of nparms*nmemb independent random normal
  # deviates, arrange into matrix of nparms, nmemb
  e <- matrix(rnorm(mean=0,sd=1,n=nmemb*nparms),
              nrow=nparms,ncol=nmemb)
  
  if(is.na(x)) {
    x <- 0
  } 
  return(t(x + D %*% e))
}

simulate_observed <- function(H,x,H_fixed=NULL,r=NULL) {
  nobs <- dim(H)[1]
  nparms <- dim(H)[2]
  if(is.null(dim(x))) {
    if(length(x)==nparms) {
      x <- matrix(x,nrow=nparms,ncol=1)
    } else {
      stop("Cannot interpret input \"x\"")
    }
  }
  nmemb <- dim(x)[2]

  retval <- H %*% x
  # Fixed components are like FF and fires, those that we do not
  # optimize
  if(!is.null(H_fixed)) {
    for(icol in 1:dim(H_fixed)[2]) {
      retval <- retval + H_fixed[,icol]
    }
  }

  if(!is.null(r)) {
    retval <- retval + rnorm(mean=0,sd=sqrt(r),n=nobs)
  }
  return(retval)
}

localization_tval <- function(dx,dy,threshold.prob=0.975) {

  nmemb <- dim(dx)[1]
  nparms <- dim(dx)[2]
  nobs <- dim(dy)[2]

  rho <- my.Pearsons.corr(t(dx),t(dy))
  tvals <- abs(rho/sqrt((1-rho^2)/(nmemb-2)))

  # Threshold Student's t value is 1.97591 for nmemb=150
  tval.thresh <- qt(p=threshold.prob,df=nmemb-1)
  retval <- matrix(1.0,nparms,nobs)

  pb <- progress.bar.start(sprintf("localizing %d obs",nobs),nobs)
  for(iobs in 1:nobs) {
    lxpar <- which(tvals[,iobs]<tval.thresh)
    if(length(lxpar)>0) {
      retval[lxpar,iobs] <- 0
    }
    pb <- progress.bar.print(pb,iobs)
  }
  progress.bar.end(pb)
  return(retval)
}

find_outliers <- function(y,dy,obs,r,rejection.threshold=3) {

  nobs <- dim(dy)[1]
  nmemb <- dim(dy)[2]

  retval <- rep(0,nobs)

  pb <- progress.bar.start(sprintf("detecting outliers for %d obs",nobs),nobs)
  
  for(iobs in 1:nobs) {
    shqhr <- sqrt((1/(nmemb-1))*sum(dy[iobs,] %*% t(dy[iobs,])) + r[iobs])
  
    # check for rejection
    z <- (y[iobs]-obs[iobs])/shqhr
    
    if(abs(z)>rejection.threshold) {
      retval[iobs] <- 1 # reject
    }
    pb <- progress.bar.print(pb,iobs)
  }
  progress.bar.end(pb)
  
  return(retval)
}

# EnSRF measurement update with localization
enkf_meas_update_loc <- function(x,dx,obs,y,dy,r,localization_mask=NULL) {

  # x   - prior state central value (nparms x 1)
  # dx  - prior state ensemble deviations (nparms x nmemb)
  # obs - observed values (nobs x 1)
  # y   - central value simulated values (nobs x 1)
  # dy  - ensemble simulated values (nobs x nmemb)
  # r   - measurement error vector (MDM); in units of variance (nobs x nobs). Not a matrix; is diagonal of formal Sz matrix
  # localization_mask - matrix of 1s and 0s to multiply K and potentially zero-out the effect of a given obs on a given parameter (nparms x nobs) 
  
  nobs <- length(y)
  nparms <- length(x)
  nmemb <- dim(dy)[2]

  sim <- y # the simulated obs from the "mean" member

    # In Whitaker & Hamill 2002 (henceforth, whitaker02a),
  # their H x'b is our dy,
  # their H Pb Ht is our HQH (Q is covariance matrix before assimilation),
  # their x'b is our dx  (parameter deviations before assimilation).

  # NOTE R's chol() returns an upper-triangular matrix, whereas
  # python's cholesky() returns a lower-triangular matrix.  We need
  # lower-triangular, so all R calls to chol() are wrapped in t().
    
  HQHR <- (1/(nmemb-1))*(t(dy) %*% dy) + diag(r)
  HQHRi <- solve(HQHR)
  
  HQ <- (1/(nmemb-1))*(t(dx) %*% dy)

  # invert central value
  
  K <- HQ %*% HQHRi
  
  # apply localization mask to K
  if(!is.null(localization_mask)) {
    K <- localization_mask * K # element-by-element multiplication
  }
  resid <- obs - sim
  x <- x + K %*% resid

  sHQHR <- t(chol(HQHR))
  sHQHRi <- solve(sHQHR)

  Ktilde <- HQ %*% t(sHQHRi) %*% solve(sHQHR + diag(sqrt(r)))  # whitaker02a, eqn 10
  
  # apply localization mask to Ktilde
  if(!is.null(localization_mask)) {
    Ktilde <- localization_mask * Ktilde # element-by-element multiplication
  }
  
  dx <- dx - t(Ktilde %*% t(dy))

  retval <- list()
  retval$dx <- dx
  retval$x <- x
  return(retval)
}

# Kalman filter measurement update.
#
# Update the solution {x,Sx} with additional 
# observations {z,Pz} obeying the form z = Hx

kf_meas_update <- function(x,Sx,H,z,Sz) {

  HPHR <- H %*% Sx %*% t(H)+Sz
  K <- (Sx %*% t(H)) %*% solve(HPHR)
  retval <- list()
  retval$x = x + K %*% (z - H %*% x)
  imkh=diag(rep(1,length(x)))-K %*% H
  
  #retval$Sx = imkh %*% Sx %*% t(imkh) + K %*% Sz %*% t(K)
  retval$Sx = imkh %*% Sx 
  
  return(retval)
}
