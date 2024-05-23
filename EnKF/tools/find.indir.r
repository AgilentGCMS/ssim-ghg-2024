find.indir <- function(indirs=c(sprintf("%s/Presentations/2024/enkf_summer_school/notebooks/SSIM-GHG",Sys.getenv("HOME")), # location on Andy's machine
                                "/work/noaa/co2/andy/Projects/enkf_summer_school/ssim-ghg-data",
                                file.path(Sys.getenv("HOME"),"shared","ssim-ghg-data")), # location on GHG Center
                       verbose=FALSE) {

  for (indir in indirs) {
    if(dir.exists(indir)) {
      if(verbose) {
        cat(sprintf("Using input dir \"%s\"\n",indir))
      }
      break
    }
  }

  if(!dir.exists(indir)) {
    stop("Error: None of the potential input directories exists.")
  }
  return(indir)
}

