# Time-stamp: <hercules-login-1.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/time-varying/persist.r: 23 May 2024 (Thu) 15:03:14 UTC>

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
indir <- find.indir()

# Load sensitivity matrices (Jacobians)
t0 <- proc.time()[3]
cat("Loading Jacobians...")
load(file.path(indir,"inversion_examples/jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
load(file.path(indir,"inversion_examples/jacobians/jacob_bgd_021624.rda"))
orig.H <- list()
orig.H$H <- jacob*(22/44) # Andrew reports units conversion needed
orig.H$H_fixed <- jacob_bgd[,c(2,3)]
rm(jacob,jacob_bgd)

cat(sprintf('%.1fs\n',proc.time()[3]-t0))

nobs <- dim(orig.H$H)[1]

# One question often asked is what happens if we get the MDM
# wrong. You can check that below, by changing r.assumed to be
# different from r.actual.
r.actual <- rep(0.5,nobs) # variance in ppm^2
#r.assumed <- rep(0.1,nobs) # variance in ppm^2
r.assumed <- r.actual

nmemb <- 3000
nparms <- 22 # 22 regions
nmons <- 24 # 24 months

Sx <- diag(rep(1,nparms))
Sx.prior <- 1e4*Sx # uninformative prior

state <- list()
x <- matrix(0,nrow=nparms,ncol=1) # prior central value
dx <- generate_ensemble(Sx=Sx.prior,nmemb=nmemb) # prior deviations

truth_condition <- matrix(NA,nrow=nmons,ncol=nparms)
for (imon in 1:nmons) {
  truth_condition[imon,] <- generate_ensemble(Sx=Sx,nmemb=1)
}

nobs <- dim(H)[1]
obs <- simulate_observed(H=H, x=truth_condition,H_fixed=H_fixed,r=r.actual)
dim(obs) <- c(nobs,1)

# Restrict to nobs randomly sampled subset of measurements. Could use
# obs_catalog or row.names of H to do more systematically-chosen
# subsets.
nobs <- 10000
lx <- sample(x=1:length(obs),size=nobs)
obs <- obs[lx,]
r.assumed <- r.assumed[lx]
r.actual <- r.actual[lx]
H <- H[lx,]
H_fixed <- H_fixed[lx,]

# generate obs and obs deviations corresponding to x and the dx
# parameter deviations.  Note that H_fixed are the fire and FF
# components sampled at the obs locations.

y <- simulate_observed(H=H,x=x,
                       H_fixed=H_fixed)
dy <- t(simulate_observed(H=H,x=t(dx),
                          H_fixed=H_fixed))

# The rejection and localization procedures are (currently) both
# serial loops over nobs and are therefore very slow. They are
# amenable to trivial parallelization, but this is not yet
# implemented. They also can be combined into one procedure, which is
# what we do operationally.

reject.outliers <- FALSE
localize <- FALSE

if(reject.outliers) {
  rej.mask <- find_outliers(y=y,dy=dy,obs=obs,r=r.assumed)

  if(any(rej.mask==1)) {
    cat(sprintf("%d observations rejected as outliers\n",
                length(which(rej.mask==1))))
    lx <- which(rej.mask==0) # these are the ones to keep
    nobs <- length(lx)
    dy <- dy[lx,]
    y <- y[lx]
    obs <- obs[lx]
    r.assumed <- r.assumed[lx]
    r.actual <- r.actual[lx]
    H <- H[lx,]
    H_fixed <- H_fixed[lx,]
    
  }
}


if(localize) {
  loc.mask <- localization_tval(dx=dx,dy=dy)
  
  cat(sprintf("%.1f%% of obs-parm relationships localized away.\n",
              100*length(which(as.vector(loc.mask)==0))/(nparms*nobs)))
} else {
  loc.mask <- NULL
}

obs_fixed <- apply(H_fixed,c(1),sum)

# obs-obs_fixed because we remove the fixed components.
enkf <- enkf_meas_update_loc(x=x,dx=dx,obs=obs-obs_fixed,
                             y=y,dy=dy,r=r.assumed,
                             localization_mask=loc.mask)

# compute sample covariance to represent posterior Sx
Sx.enkf <- cov(enkf$dx)

# Schuh invert_clean. Some different variable names,
# subset_indicator_obs can be used to subset the measurements (we do
# not use that).
subset_indicator_obs=rep(TRUE,dim(H)[1])
ret2 = invert_clean(H=H,R_diagonal=r.assumed,P_0=Sx.prior,y=obs,H_bgd=H_fixed,
                    subset_indicator_obs=subset_indicator_obs)

# Kalman filter measurement update
kf <- kf_meas_update(x=x,Sx=Sx.prior,H=H,z=obs-obs_fixed,
                     Sz=diag(r.assumed))

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

chi2.enkf <- (1/ndofs.enkf) * t(enkf$x - truth_condition) %*% solve(Sx.enkf) %*% (enkf$x - truth_condition)

chi2.ic <- (1/ndofs.ic) * t(ret2$posterior$x_hat + 1 - truth_condition) %*% solve(ret2$posterior$P) %*% (ret2$posterior$x_hat + 1 - truth_condition)

cat(sprintf("Schuh: chi2mn_state: %.2f on %d DOFs; RMSE %.2f\n",chi2.ic,ndofs.ic,compute.rmse(ret2$posterior$x_hat - truth_condition)))

cat(sprintf("   KF: chi2mn_state: %.2f on %d DOFs; RMSE %.2f\n",chi2.kf,ndofs.kf,compute.rmse(kf$x - truth_condition)))

cat(sprintf(" EnKF: chi2mn_state: %.2f on %d DOFs, RMSE %.2f\n",chi2.enkf,ndofs.enkf,compute.rmse(enkf$x - truth_condition)))

save(enkf,kf,ret2,truth_condition,file="compare_methods.results.rda")
