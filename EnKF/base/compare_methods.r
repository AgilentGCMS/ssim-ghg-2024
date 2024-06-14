# Time-stamp: <hercules-login-4.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/base/compare_methods.r: 03 Jun 2024 (Mon) 23:43:29 UTC>

# This code selects a random subset of obs, and a randomly-generated
# truth condition. Of course, the EnKF uses a randomly-generated
# ensemble. Consequently, results from a single run can be
# surprisingly different from the expected (average) performance. If
# you encounter something odd, run the script again to see it if
# persists.

source("../tools/invert_clean.r")
source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
indir <- find.indir()

# Load sensitivity matrices (Jacobians)
t0 <- proc.time()[3]
cat("Loading Jacobians...")
load(file.path(indir,"inversion_examples/jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
load(file.path(indir,"inversion_examples/jacobians/jacob_bgd_021624.rda"))
H <- jacob*(12/44) # Andrew reports units conversion needed
H_fixed <- jacob_bgd[,c(2,3)]
rm(jacob,jacob_bgd)
cat(sprintf('%.1fs\n',proc.time()[3]-t0))

nobs <- dim(H)[1]

# One question often asked is what happens if we get the MDM
# wrong. You can check that below, by changing Szd.assumed to be
# different from Szd.actual.

# "d" suffix means its the diagonal (of a diagonal matrix)
Szd.actual <- rep((0.5)^2,nobs) # variance in ppm^2
#Szd.assumed <- rep((0.1)^2,nobs) # variance in ppm^2
Szd.assumed <- Szd.actual

# For this example, EnKF approaches the best answer (that of KF) for
# an ensemble of about 8000. If the ensemble is smaller, the
# chisq-state is less than 1.0, meaning that the uncertainties are too
# large compared with the actual state errors. (Large chisq-state
# means the opposite of course--underestimated uncertainties.)
nmemb <- 8000

nparms <- 22*24 # 22 regions, 24 months

Sx <- diag(rep(1,nparms))
Sx.prior <- Sx

x.prior <- matrix(0,nrow=nparms,ncol=1) # prior central value
dx.prior <- generate_ensemble(Sx=Sx.prior,nmemb=nmemb) # prior deviations

truth_condition <- generate_ensemble(Sx=Sx,nmemb=1)
dim(truth_condition) <- c(nparms,1)

# original (non-subsetted) no. obs is defined by the number of rows in
# the H (Jacobian) matrix.
nobs <- dim(H)[1]

# Note that supplying a Szd argument to the simulate_observed function
# will result in perturbations being added to the observations.  This
# would be the place to add biases to obs.
obs <- simulate_observed(H=H, x=truth_condition,H_fixed=H_fixed,Szd=Szd.actual)
dim(obs) <- c(nobs,1)

# Restrict to nobs randomly sampled subset of measurements. Could use
# obs_catalog or row.names of H to do more systematically-chosen
# subsets.
nobs.subset <- 10000

# lx is a vector of indices into the original 1:nobs 
lx <- sample(x=1:length(obs),size=nobs.subset) 

# subset all the arrays
obs <- obs[lx,]
Szd.assumed <- Szd.assumed[lx]
Szd.actual <- Szd.actual[lx]
H <- H[lx,]
H_fixed <- H_fixed[lx,]

# We use "nobs" below, so we need to update it
nobs <- nobs.subset

# generate prior simulated obs and obs deviations corresponding to
# x.prior and dx.prior parameter deviations.  Note that H_fixed are
# the fire and FF components sampled at the obs locations.

y.prior <- simulate_observed(H=H,x=x.prior,
                             H_fixed=H_fixed)

dy.prior <- t(simulate_observed(H=H,x=t(dx.prior),
                                H_fixed=H_fixed))

# The rejection and localization procedures are (currently) both
# serial loops over nobs and are therefore very slow. They are
# amenable to trivial parallelization, but this is not yet
# implemented. They also can be combined into one procedure, which is
# what we do operationally.

reject.outliers <- FALSE
localize <- FALSE

if(reject.outliers) {
  rej.mask <- find_outliers(y=y.prior,dy=dy.prior,obs=obs,
                            r=Szd.assumed)

  if(any(rej.mask==1)) {
    cat(sprintf("%d observations rejected as outliers\n",
                length(which(rej.mask==1))))
    lx <- which(rej.mask==0) # these are the ones to keep
    nobs <- length(lx)
    dy.prior <- dy.prior[lx,]
    y.prior <- y.prior[lx]
    obs <- obs[lx]
    Szd.assumed <- Szd.assumed[lx]
    Szd.actual <- Szd.actual[lx]
    H <- H[lx,]
    H_fixed <- H_fixed[lx,]
    
  }
}


if(localize) {
  loc.mask <- localization_tval(dx=dx.prior,dy=dy.prior)
  
  cat(sprintf("%.1f%% of obs-parm relationships localized away.\n",
              100*length(which(as.vector(loc.mask)==0))/(nparms*nobs)))
} else {
  loc.mask <- NULL
}

obs_fixed <- apply(H_fixed,c(1),sum)

# obs-obs_fixed because we remove the fixed components.
enkf <- enkf_meas_update_loc(x=x.prior, dx=dx.prior,
                             obs=obs-obs_fixed, Szd=Szd.assumed,
                             y=y.prior, dy=dy.prior,
                             localization_mask=loc.mask)

# compute sample covariance to represent posterior Sx
Sx.enkf <- cov(enkf$dx)

# Schuh invert_clean. Beware of different variable names.
# subset_indicator_obs can be used to subset the measurements (we do
# not use that). His R_diagonal is assumed to be a s.d. in ppm.
subset_indicator_obs=rep(TRUE,dim(H)[1])
ic = invert_clean(H=H,R_diagonal=sqrt(Szd.assumed),
                  P_0=Sx.prior,y=obs,H_bgd=H_fixed,
                  subset_indicator_obs=subset_indicator_obs)

# Adding 1.0 to ic results because Andrew's code explicitly
# centers the scaling factors by subtracting 1.0
ic$posterior$x_hat <- ic$posterior$x_hat + 1

# Kalman filter measurement update
kf <- kf_meas_update(x=x.prior,Sx=Sx.prior,H=H,z=obs-obs_fixed,
                     Sz=diag(Szd.assumed))

posterior.dofs <- FALSE

if(posterior.dofs) {
  ndofs.ic <- ndofs.patil(ic$posterior$P)
  ndofs.kf <- ndofs.patil(kf$Sx)
  ndofs.enkf <- ndofs.patil(Sx.enkf)
} else {
  ndofs.ic <- nparms
  ndofs.kf <- nparms
  ndofs.enkf <- nparms
}

chi2.state.kf <- (1/ndofs.kf) * t(kf$x - truth_condition) %*% solve(kf$Sx) %*% (kf$x - truth_condition)
chi2.prior.kf <- (1/nparms) * t(kf$x - x.prior) %*% solve(Sx.prior) %*% (kf$x - x.prior)
obs.kf.post <- simulate_observed(H=H,x=kf$x,H_fixed=H_fixed)
chi2.obs.kf <- (1/nobs) * t(obs - obs.kf.post) %*% diag(1/Szd.assumed) %*% (obs - obs.kf.post)

chi2.state.enkf <- (1/ndofs.enkf) * t(enkf$x - truth_condition) %*% solve(Sx.enkf) %*% (enkf$x - truth_condition)
chi2.prior.enkf <- (1/nparms) * t(enkf$x - x.prior) %*% solve(Sx.prior) %*% (enkf$x - x.prior)
obs.enkf.post <- simulate_observed(H=H,x=enkf$x,H_fixed=H_fixed)
chi2.obs.enkf <- (1/nobs) * t(obs - obs.enkf.post) %*% diag(1/Szd.assumed) %*% (obs - obs.enkf.post)

chi2.state.ic <- (1/ndofs.ic) * t(ic$posterior$x_hat - truth_condition) %*% solve(ic$posterior$P) %*% (ic$posterior$x_hat - truth_condition)
chi2.prior.ic <- (1/nparms) * t(ic$posterior$x_hat - x.prior) %*% solve(Sx.prior) %*% (ic$posterior$x_hat - x.prior)
obs.ic.post <- simulate_observed(H=H,x=ic$posterior$x_hat,H_fixed=H_fixed)
chi2.obs.ic <- (1/nobs) * t(obs - obs.ic.post) %*% diag(1/Szd.assumed) %*% (obs - obs.ic.post)

#chi2.state.ic <- (1/ndofs.ic) * t(ret2$posterior$x_hat + 1 - truth_condition) %*% solve(ret2$posterior$P) %*% (ret2$posterior$x_hat + 1 - truth_condition)

cat(sprintf("   [IC] chi2 means: state %.2f, prior %.2f, obs %.2f on %d DOFs; RMSE %.2f\n",
            chi2.state.ic,chi2.prior.ic,chi2.obs.ic,ndofs.ic,compute.rmse(ic$posterior$x_hat - truth_condition)))

cat(sprintf("   [KF] chi2 means: state %.2f, prior %.2f, obs %.2f on %d DOFs; RMSE %.2f\n",
            chi2.state.kf,chi2.prior.kf,chi2.obs.kf,ndofs.kf,compute.rmse(kf$x - truth_condition)))

cat(sprintf(" [EnKF] chi2 means: state %.2f, prior %.2f, obs %.2f on %d DOFs, RMSE %.2f (%d members)\n",
            chi2.state.enkf,chi2.prior.enkf,chi2.obs.enkf,ndofs.enkf,compute.rmse(enkf$x - truth_condition),nmemb))

#save(enkf,kf,ret2,truth_condition,file="compare_methods.results.rda")
save(enkf,kf,truth_condition,file="compare_methods.results.rda")
