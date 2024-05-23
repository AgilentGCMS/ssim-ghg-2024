# Time-stamp: <hercules-login-1.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/L-curve/lcurve.r: 23 May 2024 (Thu) 16:48:26 UTC>

source("../tools/enkf.r")
source("../tools/progress.bar.r")
source("../tools/find.indir.r")
indir <- find.indir()

alphas <- 10^(-5:5)
niter <- 5

nalphas <- length(alphas)
chi2.resids <- matrix(NA,nrow=niter,ncol=nalphas)
chi2.flux <- matrix(NA,nrow=niter,ncol=nalphas)
                                      
# Load sensitivity matrices (Jacobians)
t0 <- proc.time()[3]
cat("Loading Jacobians...")
load(file.path(indir,"inversion_examples/jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
load(file.path(indir,"inversion_examples/jacobians/jacob_bgd_021624.rda"))
H <- jacob*(22/44) # Andrew reports units conversion needed
H_fixed <- jacob_bgd[,c(2,3)]
rm(jacob,jacob_bgd)
H.orig <- list()
H.orig$H <- H
H.orig$H_fixed <- H_fixed
cat(sprintf('%.1fs\n',proc.time()[3]-t0))

nobs <- dim(H)[1]

# One question often asked is what happens if we get the MDM
# wrong. You can check that below, by changing r.assumed to be
# different from r.actual.
r.actual <- rep(0.5,nobs) # variance in ppm^2
#r.assumed <- rep(0.1,nobs) # variance in ppm^2
r.assumed <- r.actual

nmemb <- 3000
nparms <- 22*24 # 22 regions, 24 months

Sx <- diag(rep(1,nparms))

n <- nalphas*niter
i <- 0

pb <- progress.bar.start(n)

ialpha <- 0
for (alpha in alphas) {
  ialpha <- ialpha + 1
  for (iter in 1:niter) {
    i <- i+1

    Sx.prior <- alpha*Sx 
    x.prior <- rep(0,nparms)
    
    truth_condition <- generate_ensemble(Sx=Sx,nmemb=1)
    dim(truth_condition) <- c(nparms,1)
    
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
    H <- H.orig$H[lx,]
    H_fixed <- H.orig$H_fixed[lx,]
    
    obs_fixed <- apply(H_fixed,c(1),sum)

    # Kalman filter measurement update
    kf <- kf_meas_update(x=x.prior,Sx=Sx.prior,H=H,z=obs-obs_fixed,
                         Sz=diag(r.assumed))

    dobs <- matrix(simulate_observed(H=H, x=kf$x,H_fixed=H_fixed) - obs,nrow=nobs,ncol=1)
    dx <- matrix(kf$x - x.prior,nrow=nparms,ncol=1)
    
    chi2.resids[iter,ialpha] <- (1/nobs) * t(dobs) %*% solve(diag(r.assumed)) %*% dobs
    chi2.flux[iter,ialpha] <- (1/nparms) * t(dx) %*% solve(Sx.prior) %*% dx

    pb <- progress.bar.print(pb,i)
    
  } # iter
} # alpha

progress.bar.end(pb)
save(chi2.resids,chi2.flux,alphas,file="lcurve.rda")
