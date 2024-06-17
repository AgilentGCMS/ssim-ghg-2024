# Time-stamp: <hercules-login-4.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/tools/enkf.r: 12 Jun 2024 (Wed) 17:05:38 UTC>

# Return data frame of the relevant pieces of a set of furnished
# obspack_id values (x)
parse.obspack_id <- function(x) {

  info <- data.frame(t(matrix(as.vector(unlist(strsplit(x,'~'))),nrow=3)),stringsAsFactors=FALSE)
  names(info) <- c("obspack_name","dataset_name","obspack_number")
  this.tempval <- data.frame(t(matrix(as.vector(unlist(strsplit(info$dataset_name,"_"))),nrow=5)),stringsAsFactors=FALSE)
  info$species <- this.tempval[,1]
  info$site <- this.tempval[,2]
  info$project <- this.tempval[,3]
  info$labcode <- this.tempval[,4]
  info$selection <- this.tempval[,5]
  return(info)
}


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

transcom.regions <- data.frame(
                               short.name=c(
                                 "NABR",
                                 "NATM",
                                 "SATR",
                                 "SATM",
                                 "NAFR",
                                 "SAFR",
                                 "EUBR",
                                 "EUTM",
                                 "ASTR",
                                 "AUST",
                                 "EURO",
                                 "NPTM",
                                 "WPTR",
                                 "EPTR",
                                 "SPTM",
                                 "NOCN",
                                 "NATM",
                                 "SATR",
                                 "SATM",
                                 "SOCN",
                                 "INTR",
                                 "INTM",
                                 "NOPT"),
                               full.name=c(
                                 "North American Boreal",
                                 "North American Temperate",
                                 "South American Tropical",
                                 "South American Temperate",
                                 "Northern Africa",
                                 "Southern Africa",
                                 "Eurasia Boreal",
                                 "Eurasia Temperate",
                                 "Tropical Asia",
                                 "Australia",
                                 "Europe",
                                 "North Pacific Temperate",
                                 "West Pacific Tropical",
                                 "East Pacific Tropical",
                                 "South Pacific Temperate",
                                 "Northern Ocean",
                                 "North Atlantic Temperate",
                                 "Atlantic Tropical",
                                 "South Atlantic Temperate",
                                 "Southern Ocean",
                                 "Indian Tropical",
                                 "South Indian Temperate",
                                 "Not optimized"),
                               type=c(rep("land",11),rep("ocean",11),"NA"),stringsAsFactors=FALSE)

# 1 T-NABR "North American Boreal"
# 2 T-NATM "North American Temperate"
# 3 T-SATR "South American Tropical"
# 4 T-SATM "South American Temperate"
# 5 T-NAFR "Northern Africa"
# 6 T-SAFR "Southern Africa"
# 7 T-EUBR "Eurasia Boreal"
# 8 T-EUTM "Eurasia Temperate"
# 9 T-ASTR "Tropical Asia"
#10 T-AUST "Australia"
#11 T-EURO "Europe"
#12 O-NPTM "North Pacific Temperate"
#13 O-WPTR "West Pacific Tropical"
#14 O-EPTR "East Pacific Tropical"
#15 O-SPTM "South Pacific Temperate"
#16 O-NOCN "Northern Ocean"
#17 O-NATM "North Atlantic Temperate"
#18 O-SATR "Atlantic Tropical"
#19 O-SATM "South Atlantic Temperate"
#20 O-SOCN "Southern Ocean"
#21 O-INTR "Indian Tropical"
#22 O-INTM "Indian Temperate"
#23 not optimized


ndofs.patil <- function(S) {
  # Compute no. of DOFs following Patil et al (2001): D. J. Patil,
  # B. R. Hunt, E. Kalnay, J. A. Yorke, and E. Ott. Local low
  # dimensionality of atmospheric dynamics. Phys. Rev. Lett.,
  # 86:5878–5881, Jun 2001. doi:10.1103/PhysRevLett.86.5878.
  #
  # Caution: not a standard measure of the DOFs
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

simulate_observed <- function(x,H,H_fixed=NULL,Szd=NULL) {
  # Note: Szd is in ppm^2 (variance), and is the diagonal of a
  # presumably-diagonal matrix.
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

  if(!is.null(Szd)) {
    retval <- retval + rnorm(mean=0,sd=sqrt(Szd),n=nobs)
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

find_outliers <- function(y,dy,obs,Szd,rejection.threshold=3) {
  # Note: Szd is in ppm^2 (variance), and is the diagonal of a
  # presumably-diagonal matrix.

  nobs <- dim(dy)[1]
  nmemb <- dim(dy)[2]

  retval <- rep(0,nobs)

  pb <- progress.bar.start(sprintf("detecting outliers for %d obs",nobs),nobs)
  
  for(iobs in 1:nobs) {
    shqhr <- sqrt((1/(nmemb-1))*sum(dy[iobs,] %*% t(dy[iobs,])) + Szd[iobs])
  
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
enkf_meas_update_loc <- function(x,dx,obs,y,dy,Szd,localization_mask=NULL) {

  # x   - prior state central value (nparms x 1)
  # dx  - prior state ensemble deviations (nparms x nmemb)
  # obs - observed values (nobs x 1)
  # y   - central value simulated values (nobs x 1)
  # dy  - ensemble deviation simulated values (nobs x nmemb)
  # Szd - measurement error vector (MDM); in units of variance (nobs x 1). Not a matrix; is diagonal of formal Sz matrix
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
    
  HQHR <- (1/(nmemb-1))*(t(dy) %*% dy) + diag(Szd)
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

  Ktilde <- HQ %*% t(sHQHRi) %*% solve(sHQHR + diag(sqrt(Szd)))  # whitaker02a, eqn 10
  
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
# observations {z,Sz} obeying the form z = Hx

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

plot.is.timeseries <- function(xs,
                               dataset_names,
                               H=NULL,
                               obs_catalog=NULL,
                               cols=NULL,
                               pdf.name=NULL) {

  # Since these load() statements can take some time, we can supply
  # them pre-loaded from memory or load them here.
  if(is.null(H)) {
    indir <- find.indir()
    t0 <- proc.time()[3]
    cat("Loading Jacobians...")
    load(file.path(indir,"inversion_examples/jacobians/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
    H <- jacob*(12/44) # Andrew reports units conversion needed
    rm(jacob)
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
  }
  
  if(is.null(obs_catalog)) {
    indir <- find.indir()
    t0 <- proc.time()[3]
    cat("Loading obs_catalog...")
    load(file.path(indir,"inversion_examples/obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda"))
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
  }

  # subset to IS only
  lx <- which(obs_catalog$TYPE=="IS")
  H <- H[lx,]
  obs_catalog <- obs_catalog[lx,]

  # change time zone to "UTC" for obs_catalog times (is by default
  # "GMT")
  attributes(obs_catalog$DATE)$tzone <- "UTC"
  
  # get dataset names from catalog
  id_components <- parse.obspack_id(obs_catalog$ID)

  if(!is.null(pdf.name)) {
    pdf(pdf.name,width=8,height=7)
    par(bg='white')
    # layout() makes subplots
    layout(matrix(1:3,nrow=3))
  } else {
    options(repr.plot.width=10,repr.plot.height=5,repr.plot.res=100)
  }

  nxs <- length(xs)
  
  # if colors not provided, generate a set
  if(is.null(cols)) {
    these.colors <- palette.colors(n = nxs, palette = "Okabe-Ito")
    cols <- list()
    for (icol in 1:nxs) {
      cols[[names(xs)[icol]]] <- these.colors[icol]
    }
  }
  
  # las=1 makes axis numbers always horizontal
  #
  # mai is the set of interior (subplot) margins in inches. Defaults
  # are too big.
  
  par(las=1,mai=c(0.5,0.6,0.6,0.1))
  
  simulated <- list()
  for (ds in dataset_names) {
    lx <- which(id_components$dataset_name==ds)
    if(length(lx)==0) {
      cat(sprintf("No such dataset \"%s\"\n",ds))
      next
    }
    nobs <- length(lx)
    cat(sprintf("%s: %d obs, from %s to %s\n",ds,nobs,
                min(obs_catalog$DATE[lx]),
                max(obs_catalog$DATE[lx])))
                
    ylim <- c(NA,NA)
    
    for (nm in names(xs)) {
      x <- xs[[nm]]

      sim <- simulate_observed(x=x,H=H[lx,])
      
      if(is.null(simulated[[ds]])) {
        simulated[[ds]] <- data.frame(nm=sim,time=obs_catalog$DATE[lx])
      } 
      simulated[[ds]][[nm]] <- sim

      ylim <- range(c(ylim,sim),na.rm=TRUE)
    }

    # cex is symbol size parameter
    cex <- 0.4
    if(nobs > 1000) {
      cex <- 0.08
    }

    # generate empty plot with predetermined axis ranges
    plot(NA,NA,
         xlim=range(simulated[[ds]]$time),
         ylim=ylim,xaxt='n',
         xlab='',ylab=expression(paste("CO"[2]," (",mu,"mol mol"^{-1},")")),
         main=ds)

    axis.POSIXct(side=1,at=seq(from=ISOdatetime(2014,9,1,0,0,0,tz="UTC"),
                               to=max(simulated[[ds]]$time),
                               by="3 months"),
                 format="%b-%Y")

    pchs <- numeric(0) # an empty numeric vector
    for (nm in c(setdiff(names(xs),"Truth"),"Truth")) { # this moves "Truth" to last in plotting order

        pch <- 20
        if(nm=="Truth") {
          cex <- 0.05
        }
      pchs <- c(pchs,pch) # append to pchs
      
      points(simulated[[ds]]$time,
           simulated[[ds]][[nm]],
           pch=pch,cex=cex,col=cols[[nm]])
    }
    
    legend(x='topright',pch=pchs,col=unlist(cols),
           legend=names(xs),horiz=TRUE,
           bty='n')
    
  }

  if(!is.null(pdf.name)) {
    dev.off() # close pdf
  }
  
  # return(sim) # an option, not necessarily recommended
}

plot.flux.timeseries <- function(ests,
                                 cols=NULL,
                                 regions=1:22,
                                 tc.priors=NULL,
                                 pdf.name=NULL) {

  # The input "ests" is a list, each element of which is a state
  # estimate (i.e., a list itself, containing x and Sx). The "truth"
  # element should only have x, not Sx.
    require(plotrix)
    
  nests <- length(ests)

  if(is.null(tc.priors)) {
    indir <- find.indir()
    t0 <- proc.time()[3]
    cat("Loading priors integrated over Transcom regions...")
    tc.priors <- load.ncdf(file.path(indir,"inversion_examples/priors/TRANSCOM_priors.nc"))
    cat(sprintf('%.1fs\n',proc.time()[3]-t0))
  }
  
  time.axis <- seq.midmon(2014,9,2016,8)
  xlim <- range(time.axis)
  
  for (nm in names(ests)) {
    # If in 1x528 space, convert to 24x22 space
    dim(ests[[nm]]$x) <- c(24,22)
    if(!is.null(ests[[nm]]$Sx)) { # because truth won't have Sx
      dim(ests[[nm]]$Sx) <- c(24*22,24*22)
    }
  }
    
  do.deviations <- FALSE
  if("Truth" %in% names(ests)) {
    do.deviations <- TRUE
    for (nm in names(ests)) {
      if(nm=="Truth") {
        next
      }
      ests[[nm]]$dx <- ests[[nm]]$x - ests$Truth$x
    }
  }

  # if colors not provided, generate a set
  if(is.null(cols)) {
    these.colors <- palette.colors(n = nests, palette = "Okabe-Ito")
    cols <- list()
    for (icol in 1:nests) {
      cols[[names(ests)[icol]]] <- these.colors[icol]
    }
  }

  if(!is.null(pdf.name)) {
    pdf(pdf.name,width=6,height=7.5)
    par(bg='white')
    # layout() makes subplots
    layout(matrix(1:4,nrow=4))
  } else {
    options(repr.plot.width=10,repr.plot.height=5,repr.plot.res=100)
  }
  
  # las=1 makes axis numbers always horizontal
  #
  # mai is the set of interior (subplot) margins in inches. Defaults
  # are too big.
  
  par(las=1,mai=c(0.5,0.6,0.6,0.1))
  cex <- 0.8
  
  for (ireg in regions) {
    
    ylim.full <- c(NA,NA)
    ylim.deviations <- c(NA,NA)
    lx.reg <- (ireg-1)*24 + 1:24
    
    for (nm in names(ests)) {

      x <- ests[[nm]]$x[,ireg]
      if(!is.null(ests[[nm]]$Sx)) {
        sx <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))
      } else {
        sx <- rep(0,24)
      }

      ylim.full <- range(c(ylim.full,(x+sx)*tc.priors$priors[ireg,],(x-sx)*tc.priors$priors[ireg,]),na.rm=TRUE)
      if(do.deviations) {
        if(nm=="Truth") {
          next
        }
        ylim.deviations <- range(c(ylim.deviations,
        	(ests[[nm]]$dx[,ireg]+sx)*tc.priors$priors[ireg,],
        	(ests[[nm]]$dx[,ireg]-sx)*tc.priors$priors[ireg,]),na.rm=TRUE)
      }
    } # nm

    plot(NA,NA,xlim=xlim,ylim=ylim.full,
         main=transcom.regions$full.name[ireg],
         xaxt='n',
         ylab='region flux (PgC/yr)',xlab='')

    axis.POSIXct(side=1,at=time.axis,format="%b-%Y")

    pchs <- numeric(0) # an empty numeric vector
    ltys <- numeric(0) # an empty numeric vector
    dts <- 30.5*86400*seq(from=-0.12,to=0.12,length.out=nests)
    iest <- 0
    for (nm in names(ests)) {
      iest <- iest+1
      dt <- dts[iest]
      flux <- ests[[nm]]$x[,ireg]*tc.priors$priors[ireg,]
      if(nm=="Truth") {
        pch <- 1
        lty <- NA
        points(x=time.axis+dt,y=flux,pch=pch,cex=cex,col=cols[[nm]])
      } else {
        pch <- 20
        lty <- 1
        flux.sd <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))*tc.priors$priors[ireg,]
        plotCI(x=time.axis+dt,y=flux,uiw=flux.sd,pch=pch,cex=cex,col=cols[[nm]],add=TRUE,sfrac=0)
      }
      pchs <- c(pchs,pch) # append to pchs
      ltys <- c(ltys,lty) # append to ltys
    }
    if(!is.null(pdf.name)) {
        y <- ylim.full[2]+0.4*diff(ylim.full)
    } else {
        y <- ylim.full[2]+0.2*diff(ylim.full)
    }

    legend(x=mean(xlim),y=y,
           pch=pchs,col=unlist(cols),xjust=0.5,
           legend=names(ests),horiz=TRUE,lty=ltys,lwd=2,
           xpd=NA,bty='n')

    if(do.deviations) {
      plot(NA,NA,xlim=xlim,ylim=ylim.deviations,
           main=sprintf("Errors for %s",transcom.regions$full.name[ireg]),
           xaxt='n',
           ylab='region flux errors (PgC/yr)',xlab='')

      axis.POSIXct(side=1,at=time.axis,format="%b-%Y")
      abline(h=0,col="grey60")

      iest <- 0
      for (nm in names(ests)) {
        iest <- iest+1
        dt <- dts[iest]
        if(nm=="Truth") {
          next
        }
        dflux <- ests[[nm]]$dx[,ireg]*tc.priors$priors[ireg,]
        dflux.sd <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))*tc.priors$priors[ireg,]
        plotCI(time.axis+dt,y=dflux,uiw=dflux.sd,pch=20,cex=cex,col=cols[[nm]],add=TRUE,sfrac=0)
      }
      these.names <- setdiff(names(cols),"Truth")
      these.cols <- cols[these.names]
      if(!is.null(pdf.name)) {
          y <- ylim.deviations[2]+0.4*diff(ylim.deviations)
      } else {
          y <- ylim.deviations[2]+0.2*diff(ylim.deviations)
      }
      legend(x=mean(xlim),
             y=y,
             pch=20,col=unlist(these.cols),xjust=0.5,
             legend=these.names,horiz=TRUE,lty=1,lwd=2,
             xpd=NA,bty='n')
    }
    
  } # ireg
  
  if(!is.null(pdf.name)) {
    dev.off() # close pdf
  }
    
}
  
 
plot.x.timeseries <- function(ests,
                              cols=NULL,
                              regions=1:22,
                              pdf.name=NULL) {

  # The input "ests" is a list, each element of which is a state
  # estimate (i.e., a list itself, containing x and Sx). The "truth"
  # element should only have x, not Sx.
  require(plotrix)
    
  nests <- length(ests)
  time.axis <- seq.midmon(2014,9,2016,8)
  xlim <- range(time.axis)
  
  for (nm in names(ests)) {
    # If in 1x528 space, convert to 24x22 space
    dim(ests[[nm]]$x) <- c(24,22)
    if(!is.null(ests[[nm]]$Sx)) { # because truth won't have Sx
      dim(ests[[nm]]$Sx) <- c(24*22,24*22)
    }
  }
    
  do.deviations <- FALSE
  if("Truth" %in% names(ests)) {
    do.deviations <- TRUE
    for (nm in names(ests)) {
      if(nm=="Truth") {
        next
      }
      ests[[nm]]$dx <- ests[[nm]]$x - ests$Truth$x
    }
  }

  # if colors not provided, generate a set
  if(is.null(cols)) {
    these.colors <- palette.colors(n = nests, palette = "Okabe-Ito")
    cols <- list()
    for (icol in 1:nests) {
      cols[[names(ests)[icol]]] <- these.colors[icol]
    }
  }

  if(!is.null(pdf.name)) {
    pdf(pdf.name,width=6,height=7.5)
    par(bg='white')
    # layout() makes subplots
    layout(matrix(1:4,nrow=4))
  } else {
    options(repr.plot.width=10,repr.plot.height=5,repr.plot.res=100)
  }
  

  # las=1 makes axis numbers always horizontal
  #
  # mai is the set of interior (subplot) margins in inches. Defaults
  # are too big.
  
  par(las=1,mai=c(0.5,0.6,0.6,0.1))
  cex <- 0.8
  
  for (ireg in regions) {
    
    ylim.full <- c(NA,NA)
    ylim.deviations <- c(NA,NA)
    lx.reg <- (ireg-1)*24 + 1:24
    
    for (nm in names(ests)) {

      x <- ests[[nm]]$x[,ireg]
      if(!is.null(ests[[nm]]$Sx)) {
        sx <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))
      } else {
        sx <- rep(0,24)
      }

      ylim.full <- range(c(ylim.full,(x+sx),(x-sx)),na.rm=TRUE)
      if(do.deviations) {
        if(nm=="Truth") {
          next
        }
        ylim.deviations <- range(c(ylim.deviations,
        	(ests[[nm]]$dx[,ireg]+sx),
        	(ests[[nm]]$dx[,ireg]-sx)),na.rm=TRUE)
      }
    } # nm

    plot(NA,NA,xlim=xlim,ylim=ylim.full,
         main=transcom.regions$full.name[ireg],
         xaxt='n',
         ylab='scaling factor',xlab='')

    axis.POSIXct(side=1,at=time.axis,format="%b-%Y")

    pchs <- numeric(0) # an empty numeric vector
    ltys <- numeric(0) # an empty numeric vector
    dts <- 30.5*86400*seq(from=-0.12,to=0.12,length.out=nests)
    iest <- 0
    for (nm in names(ests)) {
      iest <- iest+1
      dt <- dts[iest]
      flux <- ests[[nm]]$x[,ireg]
      if(nm=="Truth") {
        pch <- 1
        lty <- NA
        points(x=time.axis+dt,y=flux,pch=pch,cex=cex,col=cols[[nm]])
      } else {
        pch <- 20
        lty <- 1
        flux.sd <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))
        plotCI(x=time.axis+dt,y=flux,uiw=flux.sd,pch=pch,cex=cex,col=cols[[nm]],add=TRUE,sfrac=0)
      }
      pchs <- c(pchs,pch) # append to pchs
      ltys <- c(ltys,lty) # append to ltys
    }

      if(!is.null(pdf.name)) {
          y <- ylim.full[2]+0.4*diff(ylim.full)
      } else {
          y <- ylim.full[2]+0.2*diff(ylim.full)
      }
      
    legend(x=mean(xlim),y=y,
           pch=pchs,col=unlist(cols),xjust=0.5,
           legend=names(ests),horiz=TRUE,lty=ltys,lwd=2,
           xpd=NA,bty='n')

    if(do.deviations) {
      plot(NA,NA,xlim=xlim,ylim=ylim.deviations,
           main=sprintf("Errors for %s",transcom.regions$full.name[ireg]),
           xaxt='n',
           ylab='scaling factor errors',xlab='')

      axis.POSIXct(side=1,at=time.axis,format="%b-%Y")
      abline(h=0,col="grey60")

      iest <- 0
      for (nm in names(ests)) {
        iest <- iest+1
        dt <- dts[iest]
        if(nm=="Truth") {
          next
        }
        dflux <- ests[[nm]]$dx[,ireg]
        dflux.sd <- sqrt(diag(ests[[nm]]$Sx[lx.reg,lx.reg]))
        plotCI(time.axis+dt,y=dflux,uiw=dflux.sd,pch=20,cex=cex,col=cols[[nm]],add=TRUE,sfrac=0)
      }
      these.names <- setdiff(names(cols),"Truth")
      these.cols <- cols[these.names]
      if(!is.null(pdf.name)) {
          y <- ylim.deviations[2]+0.4*diff(ylim.deviations)
      } else {
          y <- ylim.deviations[2]+0.2*diff(ylim.deviations)
      }
      legend(x=mean(xlim),
             y=y,
             pch=20,col=unlist(these.cols),xjust=0.5,
             legend=these.names,horiz=TRUE,lty=1,lwd=2,
             xpd=NA,bty='n')
    }
    
  } # ireg
  
  if(!is.null(pdf.name)) {
    dev.off() # close pdf
  }

}
  
 
limit.z <- function(mat,zlim) {
  lx <- which(mat < zlim[1])
  if(length(lx)>1) {
    mat[lx] <- zlim[1]
  }
  lx <- which(mat > zlim[2])
  if(length(lx)>1) {
    mat[lx] <- zlim[2]
  }
  return(mat)
}
