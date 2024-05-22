# Time-stamp: <tropics.local:/local/ct/tools/R/progress.bar.r - 29 Apr 2020 (Wed) 07:22:41 MDT>

# USAGE:
#
# pb <- progress.bar.start("matching",nf)
#
# for (jf in 1:nf) { 
#
#   pb <- progress.bar.print(pb,jf)
#
#    ... stuff ...
#
# }
#
# progress.bar.end(pb)
#

progress.bar.start <- function(message,nx=NA,pct.minor=1,pct.major=10) {

  # I made a mistake when I first coded this.  Message should be
  # optional and nx required. This should look like function(nx,
  # message=NA, ... , but as you can see the arguments have message
  # first. To allow for calling this function with just the integer
  # nx, I have coded a little workaround involving a missing nx
  # argument.
  
  if(is.na(nx)) {
    nx <- as.integer(message)
    if(is.na(nx)) {
      stop("[progress.bar.start] Cannot interpret first argument as an integer.")
    }
  }
  message <- as.character(message)
  cat(sprintf("%s [",message))
  retval <- list()
  retval$start.time <- proc.time()[3]
  retval$nx <- nx
  retval$pct.minor <- pct.minor
  retval$pct.major <- pct.major
  retval$last.pct.minor.printed <- 0
  retval$last.pct.major.printed <- 0
  return(retval)
}

progress.bar.print <-  function(progress.bar,x) {

  prog.cur <- x/progress.bar$nx
  if((prog.cur - progress.bar$last.pct.major.printed)>=(progress.bar$pct.major/100)) {
    cat(sprintf('%.0f',round(prog.cur*10)))
#    cat(sprintf('[%.2f]',prog.cur*10))
    progress.bar$last.pct.major.printed <- prog.cur
    progress.bar$last.pct.minor.printed <- prog.cur
  } else {
    if((prog.cur - progress.bar$last.pct.minor.printed)>=(progress.bar$pct.minor/100)) {
      cat('.')
      progress.bar$last.pct.minor.printed <- prog.cur
    }
  }
  return(progress.bar)
     
}



new.progress.bar.print <-  function(progress.bar,x) {

  prog.cur <- x/progress.bar$nx
  cat(sprintf('%d; %.2f; %.2f; %.2f\n',x,prog.cur*10,progress.bar$last.pct.major.printed,progress.bar$pct.major/100))
  if((prog.cur - progress.bar$last.pct.major.printed)>=(progress.bar$pct.major/100)) {
    cat('o\n')
#    cat(sprintf('%d',round(prog.cur*10)))
#    cat(sprintf('[%.2f]',prog.cur*10))
    progress.bar$last.pct.major.printed <- prog.cur
    progress.bar$last.pct.minor.printed <- prog.cur
  } else {
    if((prog.cur - progress.bar$last.pct.minor.printed)>=(progress.bar$pct.minor/100)) {
      cat('.')
      progress.bar$last.pct.minor.printed <- prog.cur
    }
  }
  return(progress.bar)
     
}

progress.bar.end <- function(progress.bar=NULL,message=NULL,omit.nl=FALSE) {
  cat("] ")
  if(!is.null(progress.bar)) {
    if(!is.null(progress.bar$start.time)) {
      cat(sprintf("%.1fs",proc.time()[3]-progress.bar$start.time))
    }
  }
  if(!is.null(message)) {
    cat("  ")
    cat(message)
  }
  if(!omit.nl) {
    cat("\n")
  }
}

