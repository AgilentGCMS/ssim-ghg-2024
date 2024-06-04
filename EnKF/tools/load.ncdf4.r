# Time-stamp: <hercules-login-4.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/tools/load.ncdf4.r: 04 Jun 2024 (Tue) 00:58:46 UTC>

interpret.udunits.time <- function(vals,unitstring,tz="UTC") {
  retval <- list()
  retval$is.time <- FALSE
  if(length(grep("^DAYS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    # length of string is a bad parsing heuristic, but works
    # for everything so far.
    #     21 chars if "days since 1900-01-01"
    #     19 chars if "days since 1900-1-1"
    retval$vals <- as.POSIXct(substr(unitstring,11,nchar(unitstring)),
                              tz=tz) + vals*86400
    retval$tz <- tz # UTC tzone is a presumption
  }
  # ECCO-Darwin time units are "number of days from January 0, 0000", but units are "MATLAB datenum"
  if(length(grep("^MATLAB datenum$", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    t0 <- ISOdatetime(0,1,1,0,0,0,tz="UTC")-86400  
    retval$vals <- t0 + vals*86400
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^SECONDS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- as.POSIXct(substr(unitstring,15,nchar(unitstring)),
                              tz=tz) + vals
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^HOURS since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- as.POSIXct(substr(unitstring,12,nchar(unitstring)),
                              tz=tz) + vals*3600
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^MINUTES since", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- as.POSIXct(substr(unitstring,15,nchar(unitstring)),
                              tz=tz) + vals*60
    retval$tz <- tz # UTC tzone is a presumption
  }
  if(length(grep("^decimal date", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- decimal.to.POSIX(vals)
    retval$tz <- tz # UTC tzone is a presumption
  }
  return(retval)
}
  
load.ncdf <- function(ncname,lowercase=FALSE,dims=TRUE,attribs=NULL,vars=NULL,
                      verbose=FALSE,quiet=FALSE,
                      collapse_degen=TRUE) {    # collapse_degen new, 12 Aug 2016 ARJ

  # if attribs is NULL, load no global attributes
  # if attribs is NA, load all global attributes
  
  if(verbose) {
    cat(sprintf("[load.ncdf]  Called with ncname = %s,\n",ncname))
    cat(sprintf("[load.ncdf]  and lowercase = %s,\n",lowercase))
    cat(sprintf("[load.ncdf]  and dims = %s,\n",dims))
    cat(sprintf("[load.ncdf]  and collapse_degen = %s,\n",collapse_degen))
    if(!is.null(attribs)) {
      cat(sprintf("[load.ncdf]  and attribs = %s,\n",attribs))
    }
    if(!is.null(vars)) {
      cat(sprintf("[load.ncdf]  and vars = %s.\n",vars))
    }
  }
  
  
  library(ncdf4)
  
  retval <- list()

  if(!verbose) {
    sink(file="/dev/null")
  }
  nc <- nc_open(ncname)
  i <- 0

  # variables
  if(is.null(vars)) {
    vars <- names(nc$var)
  }


  for(var in vars) {
    if(! var %in% names(nc$var)) {
      if(!quiet) {
        warning(sprintf("variable \"%s\" not in file \"%s\".",var,ncname))
      }
      next
    }
    i <- i+1
    if(!quiet) {
      cat(sprintf("loading variable %d, \"%s\"\n",i,var))
    }
    if(nc$var[[var]]$prec=="char") {
      # 3 Nov 3026, experiencing trouble with character strings and collapse_degen=FALSE
      retval[[i]] <- ncvar_get(nc,var,collapse_degen=TRUE)  
    } else {
      retval[[i]] <- ncvar_get(nc,var,collapse_degen=collapse_degen)  # collapse_degen new, 12 Aug 2016 ARJ
    }
    
    # retrieve all attributes
    vatts <- ncatt_get(nc,var)
    attributes(retval[[i]]) <- c(attributes(retval[[i]]),vatts)

    if("units" %in% names(vatts)) {
      time.checked <- interpret.udunits.time(vals=retval[[i]],unitstring=vatts$units)
      if(time.checked$is.time) {
        retval[[i]] <- time.checked$vals
        attributes(retval[[i]])$tzone <- time.checked$tz 
      }
    }
    
    if(FALSE) {
      # retrieve units, if available
      vunits <- ncatt_get(nc,var,attname="units")
      if(vunits$hasatt) {
        attributes(retval[[i]])$units <- vunits$value
        time.checked <- interpret.udunits.time(vals=retval[[i]],unitstring=vunits$value)
        if(time.checked$is.time) {
          retval[[i]] <- time.checked$vals
          attributes(retval[[i]])$tzone <- time.checked$tz 
        }
        
      }
      
      # retrieve long_name, if available
      vunits <- ncatt_get(nc,var,attname="long_name")
      if(vunits$hasatt) {
        attributes(retval[[i]])$long_name <- vunits$value
      }
    }
# Aug 2011, ARJ:  _FillValue is automatically handled by ncdf4 package.
    
#    # check for _FillValue (as opposed to missval)
#    fv <- ncatt_get(nc,var,attname="_FillValue")
#    if(fv$hasatt) {
#      lx <- which(retval[[i]] == fv$value)
#      if(verbose) {
#        cat(sprintf('[load.ncdf] _FillValue detected for %d values in variable %s.\n',
#                    length(lx),var))
#      }
#      retval[[i]][lx] <- NA
#    }
#    fv <- ncatt_get(nc,var,attname="_fillValue")
#    if(fv$hasatt) {
#      lx <- which(retval[[i]] == fv$value)
#      if(verbose) {
#        cat(sprintf('[load.ncdf] _fillValue detected for %d values in variable %s.\n',
#                    length(lx),var))
#      }
#      retval[[i]][lx] <- NA
#    }
    # 10 Sep 2009:  scalars have ndims==0, so the following is broken.
#    if(nc$var[[var]]$ndims > 0) {
#      retval[[i]] <- get.var.ncdf(nc,var)
#    } else {  # empty variable (present in nc files from ORNL)
#      retval[[i]] <- NA
#    }
    if(lowercase) {
      names(retval)[i] <- gsub('_','.',var)
    } else {
      names(retval)[i] <- var
    }
  }

  # dimensions
  if(dims) {
    retval$dim <- nc$dim
    for (d in names(nc$dim)) {
      ddot <- d
      if(lowercase) {
        ddot <- gsub('_','.',ddot)
      }
      retval[[ddot]] <- nc$dim[[d]]$vals
      unitstring <- nc$dim[[d]]$units

      # retrieve all attributes
      vatts <- ncatt_get(nc,d)
      attributes(retval[[ddot]]) <- c(attributes(retval[[ddot]]),vatts)
      
      time.checked <- interpret.udunits.time(nc$dim[[d]]$vals,nc$dim[[d]]$units)
      if(time.checked$is.time) {
        retval[[ddot]] <- time.checked$vals
        attributes(retval[[ddot]])$tzone <- time.checked$tz # UTC tzone is a presumption
      }
      

    }
  }

  if(!is.null(attribs)) {
    if((length(attribs)==1) && (is.na(attribs))) {
      attribs <- names(ncatt_get(nc,0,NA))
    }
  }
  for (att in attribs) {
    target.name <- att
    if(lowercase) {
      target.name <- gsub('_','.',target.name)
    }
    # Some attribute names are reserved and we should
    # not overwrite them with values from the netCDF
    # file.  An example is "names", which occurs in
    # column3d output.  Tack on a "2" to such names.
    if (target.name %in% names(attributes(retval))) {
      target.name <- paste(target.name,"2",sep='')
    }
    
    attributes(retval)[[target.name]] <- NULL
    att.val <- ncatt_get(nc,0,att)
    if(att.val$hasatt) {
      attributes(retval)[[target.name]] <- att.val$value
    }  else {
      if(!quiet) {
        warning(sprintf("[load.ncdf] Cannot find attribute \"%s\" in file \"%s\".",att,nc$filename))
      }
    }
  }

  nc_close(nc)
  if(!verbose) {
    sink(NULL)
  }
  if(lowercase) {
    names(retval) <- tolower(names(retval))
  }
  return(retval)
}
