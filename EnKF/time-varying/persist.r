# Time-stamp: <hercules-login-1.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/toy/toy.r: 22 May 2024 (Wed) 20:19:26 UTC>

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
indir <- find.indir()

###############################################
#--  Load sensitivity matrices
###############################################
load(file.path(indir,"data/sensitivity_matrices/","trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
load(file.path(indir,"data/sensitivity_matrices/","jacob_bgd_021624.rda"))


R_diagonal_in = rep(.01,(dim(jacob)[1])) # variance in ppm^2


nmemb <- 5000
nparms <- 22 # 22 regions
nmon <- 24

Sx <- diag(rep(1,nparms))
Sx.prior <- 1e4*Sx # uninformative prior

x <- matrix(0,nrow=nparms,ncol=1) # prior central value
dx <- generate_ensemble(Sx=Sx.prior,nmemb=nmemb) # prior deviations

truth_condition <- matrix(NA,nrow=nmon,ncol=nparms)
for (imon in 1:nmon) {
  truth_condition[imon,] <- generate_ensemble(Sx=Sx,nmemb=1)
}

nobs <- dim(jacob)[1]
r <- rep(3,nobs)
obs <- simulate_observed(H=jacob, x=truth_condition,H_fixed=jacob_bgd[,c(2,3)],r=r)
dim(obs) <- c(nobs,1)

# number of obs per month
nobs <- 100

# retrict to nobs subset
lx <- sample(x=1:length(obs),size=nobs)
#dy <- dy[lx,]
#y <- y[lx]
obs <- obs[lx,]
r <- r[lx]
jacob <- jacob[lx,]
jacob_bgd <- jacob_bgd[lx,]


# generate obs and obs deviations corresponding to x and the dx
# parameter deviations.  Note that H_fixed are the fire and FF
# components sampled at the obs locations.

y <- simulate_observed(H=jacob,x=x,
                       H_fixed=jacob_bgd[,c(2,3)])
dy <- t(simulate_observed(H=jacob,x=t(dx),
                          H_fixed=jacob_bgd[,c(2,3)]))



# The rejection and localization procedures are (currently) both
# serial loops over nobs and are therefore very slow. They are
# amenable to trivial parallelization, but this is not yet
# implemented. They also can be combined into one procedure, which is
# what we do operationally.

reject.outliers <- FALSE
localize <- FALSE

if(reject.outliers) {
  rej.mask <- find_outliers(y=y,dy=dy,obs=obs,r=r)

  if(any(rej.mask==1)) {
    cat(sprintf("%d observations rejected as outliers\n",
                length(which(rej.mask==1))))
    lx <- which(rej.mask==0) # these are the ones to keep
    dy <- dy[lx,]
    y <- y[lx]
    obs <- obs[lx]
    r <- r[lx]
    H <- jacob[lx,]
    H_fixed <- jacob_bgd[lx,]
  }
}

nobs <- length(obs)

if(localize) {
  loc.mask <- localization_tval(dx=dx,dy=dy)
  
  cat(sprintf("%.1f%% of obs-parm relationships localized away.\n",
              100*length(which(as.vector(loc.mask)==0))/(nparms*nobs)))
} else {
  loc.mask <- NULL
}
obs_bgd <- apply(jacob_bgd[,c(2,3)],c(1),sum)

# obs-obs_bgd because we remove the fixed components.
post <- enkf_meas_update_loc(x=x,dx=dx,obs=obs-obs_bgd,
                             y=y,dy=dy,r=r,
                             localization_mask=loc.mask)

Sx.enkf <- cov(post$dx)

subset_indicator_obs=rep(TRUE,dim(jacob)[1])
ret2 = invert_clean(H=jacob,R_diagonal=r,P_0=Sx.prior,y=obs,H_bgd=jacob_bgd,
                    subset_indicator_obs=subset_indicator_obs)

kf <- kf_meas_update(x=x,Sx=Sx.prior,H=jacob,z=obs-obs_bgd,
                     Sz=diag(r))

posterior.dofs <- TRUE

if(posterior.dofs) {
  ndofs.ic <- ndofs.patil(ret2$posterior$P)
  ndofs.kf <- ndofs.patil(kf$Sx)
  ndofs.enkf <- ndofs.patil(Sx.enkf)
} else {
  ndofs.ic <- nparms
  ndofs.kf <- nparms
  ndofs.enkf <- nparms
}

chi2.kf <- (1/ndofs.kf) * t(kf$x - truth_condition) %*% solve(kf$Sx) %*% (kf$x - truth_condition)

chi2.enkf <- (1/ndofs.enkf) * t(post$x - truth_condition) %*% solve(Sx.enkf) %*% (post$x - truth_condition)

chi2.ic <- (1/ndofs.ic) * t(ret2$posterior$x_hat + 1 - truth_condition) %*% solve(ret2$posterior$P) %*% (ret2$posterior$x_hat + 1 - truth_condition)

cat(sprintf("Schuh: chi2mn_state: %.2f on %d DOFs; RMSE %.2f\n",chi2.ic,ndofs.ic,compute.rmse(ret2$posterior$x_hat - truth_condition)))

cat(sprintf("   KF: chi2mn_state: %.2f on %d DOFs; RMSE %.2f\n",chi2.kf,ndofs.kf,compute.rmse(kf$x - truth_condition)))

cat(sprintf(" EnKF: chi2mn_state: %.2f on %d DOFs, RMSE %.2f\n",chi2.enkf,ndofs.enkf,compute.rmse(post$x - truth_condition)))

save(post,kf,ret2,truth_condition,file="res.rda")
