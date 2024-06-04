# Time-stamp: <hercules-login-4.hpc.msstate.edu:/work/noaa/co2/andy/Projects/enkf_summer_school/repo/ssim-ghg-2024/EnKF/tools/time.r: 03 Jun 2024 (Mon) 18:10:35 UTC>

idate.to.POSIX <- function(idate,tz="UTC") {
  # idate is a TM5 construct of 6 integers representing a UTC
  # date-time. This construct cannot consider fractional seconds.
  # takes a vector of integers representing ONE date cannot now handle
  # multiple dates year month day required (00:00:00 assumed) if idate
  # has 6 elements, will use 4,5,6 as HH::MM::SS UTC
  if(length(idate)==3) {
    return(ISOdatetime(idate[1],idate[2],idate[3],0,0,0,tz=tz))
  }
  if(length(idate)==6) {
    return(ISOdatetime(idate[1],idate[2],idate[3],idate[4],idate[5],idate[6],tz=tz))
  }
  print(idate)
  stop(sprintf("Was expecting 3 or 6 elements in the idate numeric vector, got %d instead.",length(idate)))
  
}

POSIX.to.idate <- function(dates,tz="UTC") {
  lt <- as.POSIXlt(dates)
  return(data.frame(year=lt$year+1900,
                    month=lt$mon+1,
                    day=lt$mday,
                    hour=lt$hour,
                    minute=lt$min,
                    second=lt$sec))
}

itau.to.POSIX <- function(itaus,year0=1999) {
  # itau is a TM5 construct. Like POSIX times, it is
  # seconds-since-epoch, but the epoch tends to be different for
  # different TM5 runs. Recall of course that the UNIX epoch is
  # 1970-01-01 00:00:00 Z
  if(is.na(year0)) {
    year0 <- 1999
  }
  return(ISOdatetime(year0,1,1,0,0,0,tz="UTC")+itaus)
}

itau2date <- itau.to.POSIX

is.leap <- function(year) {
    return(year%%4 == 0 & (year%%100 != 0 | year%%400 == 0))
}

days.in.year <- function(years) {
  retval <- numeric(0)
  for (yr in years) {
    if(is.leap(yr)) {
      retval <- c(retval,366)
    } else {
      retval <- c(retval,365)
    }
  }
  return(retval)
}

days.in.month <- function(year=NULL) {
  # Returns a vector of 12 days-per-month for the specified year.  If
  # year is not specified, just return non-leap years days per month.
  retval <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  if(!is.null(year)) {
    if(is.leap(year)) {
      retval[2] <- 29
    }
  }
  return(retval)
}


seconds.in.year <- function(year) {
    # Returns the integer number of seconds in the specified
    # year. Leap seconds are not considered.  See chron package for
    # better handling.
    return(86400*days.in.year(year))
}

# (my.)length.POSIXlt is a candidate for deletion since this
# capability is now provided by base.
my.length.POSIXlt <- function(x) {
  # It's pointlessly difficult to retrieve the number of dates
  # in a POSIXlt vector.  This is how you do it.
  return(length(x$wday))
}

format.et <- function(elapsed.seconds,verbose=FALSE) {
  # Convenient way to provide a compact character string summarizing
  # the elapsed time.
  tvec <- character(0)
  years <- 0
  weeks <- 0
  days <- 0
  hours <- 0
  minutes <- 0
  seconds <- 0
  
  # 365-day years
  while(elapsed.seconds > 86400*365) {
    elapsed.seconds <- elapsed.seconds - 86400*365
    years <- years + 1
  }

  # weeks
  while(elapsed.seconds > 86400*7) {
    elapsed.seconds <- elapsed.seconds - 86400*7
    weeks <- weeks + 1
  }

  # days
  while(elapsed.seconds > 86400) {
    elapsed.seconds <- elapsed.seconds - 86400
    days <- days + 1
  }

  # hours
  while(elapsed.seconds > 3600) {
    elapsed.seconds <- elapsed.seconds - 3600
    hours <- hours + 1
  }

  # minutes
  while(elapsed.seconds > 60) {
    elapsed.seconds <- elapsed.seconds - 60
    minutes <- minutes + 1
  }

  # seconds <- elapsed.seconds
  if(verbose) {
    if(years>0) {
      tvec <- c(tvec,sprintf("%d year%s",years,ifelse(years>1,"s","")))
    }
    if(weeks>0) {
      tvec <- c(tvec,sprintf("%d week%s",weeks,ifelse(weeks>1,"s","")))
    }
    if(days>0) {
      tvec <- c(tvec,sprintf("%d day%s",days,ifelse(days>1,"s","")))
    }
    if(hours>0) {
      tvec <- c(tvec,sprintf("%d hour%s",hours,ifelse(hours>1,"s","")))
    }
    if(minutes>0) {
      tvec <- c(tvec,sprintf("%d minute%s",minutes,ifelse(minutes>1,"s","")))
    }
    if(elapsed.seconds>0) {
      tvec <- c(tvec,sprintf("%.0f second%s",elapsed.seconds,ifelse(elapsed.seconds>1,"s","")))
    }
  } else {
    if(years>0) {
      tvec <- c(tvec,sprintf("%dy",years))
    }
    if(weeks>0) {
      tvec <- c(tvec,sprintf("%dw",weeks))
    }
    if(days>0) {
      tvec <- c(tvec,sprintf("%dd",days))
    }
    if(hours>0) {
      tvec <- c(tvec,sprintf("%02dh",hours))
    }
    if(minutes>0) {
      tvec <- c(tvec,sprintf("%02dm",minutes))
    }
    if(elapsed.seconds>0) {
      tvec <- c(tvec,sprintf("%02.0fs",elapsed.seconds))
    }
  }
  return(paste(tvec,collapse='-'))
}


decimal.to.POSIX <- function(decimal.dates,tz="UTC") {

  # converts a vector of decimal dates (floating point numbers) into
  # POSIXct objects

  retval <- NULL
  
  for (dd in decimal.dates) {
    year <- trunc(dd)
    y0 <- ISOdatetime(year,1,1,0,0,0,tz="UTC")
    y1 <- as.numeric(ISOdatetime(year+1,1,1,0,0,0,tz="UTC"))
    sec.in.year <- y1 - as.numeric(y0)
    pd <- as.POSIXct(y0 + (dd-year)*sec.in.year,tz=tz)
    if(is.null(retval)) {
      retval <- pd
    } else {
      retval <- c(retval,pd)
    }
  }
  attributes(retval)$tzone <- tz
  
  return(retval)
}


POSIX.to.decimal <- function(posix.dates,tz="UTC") {

  # convert a vector of POSIXct dates into floating point decimal
  # dates.
  
  if(!any(class(posix.dates)=="POSIXt")) {
    stop("[decimal.year] Expecting one of class(posix.dates) == 'POSIXt'")
  }

  # Code as of 3 Aug 2009: much faster than previous implementation,
  # as there is no looping.

  #   There's a potential problem here in that I'm not sure how the R
  #   date-time classes handle a vector of times with different time
  #   zone attributes among the members.  It would appear that the
  #   attribute belonging to the vector overrides any such attribute
  #   assigned to a member.  This is probably because it is the vector
  #   object which can have attributes, not the member.  This code
  #   assumes that all members of the vector have the same time zone
  #   attribute.

  # Note that the POSIXlt "year" is actual year - 1900
  lt <- as.POSIXlt(posix.dates,tz=tz) # vector of POSIX local time
  
  # y0 and y1 are vectors, in seconds, of the start of year n and n+1.
  y0 <- as.numeric(ISOdatetime(lt$year+1900,1,1,0,0,0,tz=attributes(lt)$tzone))
  y1 <- as.numeric(ISOdatetime(lt$year+1901,1,1,0,0,0,tz=attributes(lt)$tzone))
  return(lt$year+1900+(as.numeric(posix.dates)-y0)/(y1-y0))

}

# POSIXct stores time internally as the signed number of seconds since
# the Unix epoch 1970-01-01 00:00:00 UTC.  See help(POSIXt).  These
# epoch.seconds functions convert between this numeric value and the
# POSIXct object.

epoch.seconds.to.POSIX <- function(sec,tzone="UTC") {
  class(sec) <- c("POSIXt","POSIXct")
  attributes(sec)$tzone <- tzone
  return(sec)
}

POSIX.to.epoch.seconds <- function(pt) {
  return(as.numeric(pt))
}

seq.midmon <- function(year0,month0,year1,month1) {

  # Return a vector of UTC mid-month POSIXct values, spaced unequally
  # as they should be.  This is useful to make a time axis for monthly
  # data. N.B. seq.POSIXt doesn't do this right.

  retval <- NULL
  y <- year0
  m <- month0

  while(TRUE) {
    if(is.null(retval)) {
      retval <- month.middle(y,m)
    } else {
      retval <- c(retval,month.middle(y,m))
    }
    if(y == year1 & m == month1) {
      break
    }
    m <- m+1
    if(m>12) {
      y <- y+1
      m <- 1
    }
  }
  attributes(retval)$tzone <- "UTC"
  return(retval)
}

month.middle <- function(year,month) {
  dpm <- days.in.month(year)[month]  # can have leap years
  if((year == 1969) & (month==12)) {
    retval <- (mean(c(ISOdatetime(year,month,1,0,0,0,tz="UTC"),
                      epoch.seconds.to.POSIX(-1))))
  } else {
    retval <- (mean(c(ISOdatetime(year,month,1,0,0,0,tz="UTC"),
                      ISOdatetime(year,month,dpm,23,59,59,tz="UTC"))))
  }
  attributes(retval)$tzone <- "UTC"
  return(retval)
}

seconds.in.month <- function(year,month) {
  return(86400*days.in.month(year)[month])  # can have leap years
}

#Julian Date in decimal notation.  

julian.to.POSIX <- function(julians,epoch=ISOdatetime(0001,1,1,0,0,0,tz="UTC")) {
  # Convention is that the Julian date of the epoch is 1.0 (not 0!)
  retval <- epoch+86400*(julians-1.0)
  attributes(retval)$tzone <- attributes(epoch)$tzone
  return(retval)
}

POSIX.to.julian <- function(times,epoch=ISOdatetime(0001,1,1,0,0,0,tz="UTC")) {
  # Convention is that the Julian date of the epoch is 1.0 (not 0!)
  return(as.numeric(difftime(times,epoch,units="days")+1.0))
}


