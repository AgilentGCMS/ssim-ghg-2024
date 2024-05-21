#-- Create the true flux field from prior and true state vector

pull_true_transcom_flux = function(prior_flux_file,state_true)
{
  con = nc_open(prior_flux_file)
  NEE = ncvar_get(con,"NEE")
  nc_close(con)

  NEE_1x1 = aaply(NEE,3,.fun=function(x){expand_5x4_2_1x1(x)}) %>%
           aperm(c(2,3,1))

  zz = aaply(NEE_1x1,3,.fun=function(x){grid2transcom(x,file_location=data_dir)})
  yy = matrix(state_true,nrow=24,byrow=FALSE)
  qq = zz[,-1]*yy
  transcom_fluxes_real = apply(qq,2,sum) *30.5*3600*24*1e3*0.5 * 1e-15 # ~ PgC/yr adjustment to prior
  ret = list(monthly=qq,annual_2yr=transcom_fluxes_real)
  return(ret)
}

generate_observations = function(H,H_bgd,state_vector,err_obs=NULL){
  obs = H %*% (1+state_vector) #+ apply(H_bgd[,c(2,3)],1,sum)
  if(!is.null(err_obs)){obs = obs + rnorm(err_obs)}
  return(obs)
}


grid2transcom = function(mat,model.grid.x=1,model.grid.y=1,transcom=TRUE,
                         file_location=data_dir)
{
  tfile = paste(file_location,"transcom/TRANSCOM_mask_GEOS_Chem_",model.grid.y,"x",model.grid.x,".nc",sep="")
  afile = paste(file_location,"areas/area_",model.grid.x,"x",model.grid.y,".nc",sep="")

  grid_areas = load.ncdf(afile,vars=c("grid_cell_area"))$grid.cell.area
  
  transfil = load.ncdf(tfile)
  if(model.grid.x==1 & model.grid.y==1){
    regions = transfil$mask64
  }else{
      for(i in 1:23){
        if(i==1){ret = transfil[[1]]}else{
          transfil[[i]][transfil[[i]]==1] = i
          ret = ret + transfil[[i]]}
      }
  regions=ret
  }

  #-- weird 3.63216561952993e-40 value in transcom files
  regions[regions > 0 & regions < 1e-6] = 0
  
 if(transcom)
  {
    transcom.agg  = aggregate(as.vector(mat * grid_areas),list(as.vector(regions)),sum)
    transcom.agg.vect = transcom.agg[,2] #as.vector(unlist(transcom.agg))
    
    return(transcom.agg.vect)
  }else{
    return(mat * grid_areas)
  }
}



transcom2grid = function(vect=rep(2,23),model.grid.x=1,model.grid.y=1,transcom=TRUE,
                         file_location=data_dir)
{
  tfile = paste(file_location,"transcom/TRANSCOM_mask_GEOS_Chem_",model.grid.y,"x",model.grid.x,".nc",sep="")
  afile = paste(file_location,"areas/area_",model.grid.x,"x",model.grid.y,".nc",sep="")
  
  if(!file.exists(tfile)){stop(paste("file:",tfile," doesn't exist"))}
  if(!file.exists(afile)){stop(paste("file:",afile," doesn't exist"))}
  
  grid_areas = load.ncdf(afile,vars=c("grid_cell_area"))$grid.cell.area

  transfil = load.ncdf(tfile)
  if(model.grid.x==1 & model.grid.y==1){
    regions = transfil$mask64
  }else{
    for(i in 1:23){
      if(i==1){ret = transfil[[1]]}else{
        transfil[[i]][transfil[[i]]==1] = i
        ret = ret + transfil[[i]]}
    }
    regions=ret
  }
  
  regions[regions < 0.001] = 0
  
  for(i in 1:22)
  {
    regions[regions == i] = vect[i]
  }
  
 return(regions)
}

#-- checking effect of truncating big and small scalings
tm = function(x,low,high){
  x[x < low] = low;x[x > high] = high
  return(x)
}

expand_5x4_2_1x1 = function(x){
  v = as.vector(x)
  m = matrix(as.vector(matrix(rep(v,5),nrow=5,byrow=TRUE)),nrow=360)
  m2 = t(m)
  v2 = as.vector(m2)
  m3 = matrix(as.vector(matrix(rep(v2,4),nrow=4,byrow=TRUE)),nrow=184)
  m4 = t(m3)
  return(m4[,3:182])
}

#-- EXTRA CODE
#gridded_5x4_state = transcom2grid(as.vector(inv_object$x_hat),model.grid.x = 5,model.grid.y = 4,
#                                  file_location = tr_dir)
#
#gridded_5x4_mean_flux = aaply(NEE,c(3),.fun=function(x){x*gridded_5x4_state}) %>%
#  aperm(c(2,3,1))
#
#gridded_5x4_state_samples = aaply(samps,1,.fun=function(x){transcom2grid(x,model.grid.x = 5,model.grid.y = 4,file_location = tr_dir)}) %>%
#  aperm(c(2,3,1))
#
#gridded_5x4_flux_samples = array(NA,dim=c(dim(gridded_5x4_mean_flux),sample_number))
#
#for(i in 1:sample_number){
#  tmp_state = gridded_5x4_state_samples[,,i]
#  gridded_5x4_flux_samples[,,,i] = aaply(gridded_5x4_mean_flux,3,.fun=function(x){x*tmp_state}) %>%
#    aperm(c(2,3,1))
#}

# Time-stamp: <andyj.cmdl.noaa.gov:/ct/tools/R/load.ncdf4.r - 07 Nov 2011 (Mon) 15:26:41 MST>

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
  if(length(grep("^decimal date", unitstring,ignore.case=TRUE))==1) {
    retval$is.time <- TRUE
    retval$vals <- decimal.to.POSIX(vals)
    retval$tz <- tz # UTC tzone is a presumption
  }
  return(retval)
}


load.ncdf <- function (ncname, lowercase = TRUE, dims = TRUE, attribs = NULL, 
                       vars = NULL, verbose = FALSE) 
{
  if (verbose) {
    cat(sprintf("[load.ncdf]  Called with ncname = %s,\n", 
                ncname))
    cat(sprintf("[load.ncdf]  and lowercase = %s,\n", lowercase))
    cat(sprintf("[load.ncdf]  and dims = %s,\n", dims))
    if (!is.null(attribs)) {
      cat(sprintf("[load.ncdf]  and attribs = %s,\n", attribs))
    }
    if (!is.null(vars)) {
      cat(sprintf("[load.ncdf]  and vars = %s.\n", vars))
    }
  }
  library(ncdf4)
  retval <- list()
  nc <- nc_open(ncname)
  i <- 0
  if (is.null(vars)) {
    vars <- names(nc$var)
  }
  for (var in vars) {
    if (!var %in% names(nc$var)) {
      warning(sprintf("variable \"%s\" not in file \"%s\".", 
                      var, ncname))
      next
    }
    i <- i + 1
    retval[[i]] <- ncvar_get(nc, var)
    vunits <- ncatt_get(nc, var, attname = "units")
    if (vunits$hasatt) {
      attributes(retval[[i]])$units <- vunits$value
      time.checked <- interpret.udunits.time(vals = retval[[i]], 
                                             unitstring = vunits$value)
      if (time.checked$is.time) {
        retval[[i]] <- time.checked$vals
        attributes(retval[[i]])$tzone <- time.checked$tz
      }
    }
    if (lowercase) {
      names(retval)[i] <- gsub("_", ".", var)
    }
    else {
      names(retval)[i] <- var
    }
  }
  if (dims) {
    retval$dim <- nc$dim
    for (d in names(nc$dim)) {
      ddot <- d
      if (lowercase) {
        ddot <- gsub("_", ".", ddot)
      }
      retval[[ddot]] <- nc$dim[[d]]$vals
      unitstring <- nc$dim[[d]]$units
      time.checked <- interpret.udunits.time(nc$dim[[d]]$vals, 
                                             nc$dim[[d]]$units)
      if (time.checked$is.time) {
        retval[[ddot]] <- time.checked$vals
        attributes(retval[[ddot]])$tzone <- time.checked$tz
      }
    }
  }
  for (att in attribs) {
    target.name <- att
    if (lowercase) {
      target.name <- gsub("_", ".", target.name)
    }
    if (target.name %in% names(attributes(retval))) {
      target.name <- paste(target.name, "2", sep = "")
    }
    attributes(retval)[[target.name]] <- NULL
    att.val <- ncatt_get(nc, 0, att)
    if (att.val$hasatt) {
      attributes(retval)[[target.name]] <- att.val$value
    }
    else {
      warning(sprintf("[load.ncdf] Cannot find attribute \"%s\" in file \"%s\".", 
                      att, nc$filename))
    }
  }
  nc_close(nc)
  if (lowercase) {
    names(retval) <- tolower(names(retval))
  }
  return(retval)
}

my.col = function (nl) 
{
  rgb(c(seq(0, 1, length = floor(nl/2)), rep(1, ceiling(nl/2))), 
      c(abs(sin((0:(nl - 1)) * pi/(nl - 1)))), c(c(rep(1, ceiling(nl/2)), 
                                                   seq(1, 0, length = floor(nl/2)))))
}

#-- Plotting of monthly flux averages for 9/2014 - 8/2016 for chosen Transcom Regions
#- need to make general so next line can take vector of region #s

plot_transcom_flux_by_month = function(data=ret){
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #transcom_region_plot_monthly = 1
  #$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  #-- Plot transcom regions 
  for(transcom_region_plot_monthly in 1:22){
    subset = ret$prior_tr_monthly[,,transcom_region_plot_monthly]
    prior_df = cbind(FLUX=as.vector(subset),KIND=rep("Prior",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
    
    subset = ret$post_tr_monthly[,,transcom_region_plot_monthly]
    post_df = cbind(FLUX=as.vector(subset),KIND=rep("Post",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
    
    combined_df  = as.data.frame(rbind(post_df,prior_df))
    
    combined_df$KIND = factor(combined_df$KIND)
    combined_df$FLUX = as.numeric(combined_df$FLUX)
    combined_df$MONTH = factor(combined_df$MONTH,levels=1:24)
    require(ggplot2)
    
    options(repr.plot.width=20, repr.plot.height=8)
    
    #-- This plot is monthly avg regional flux (in units of PgC/yr), plotted for 2 years
    
    subset_truth = ret$transcom_fluxes_real_monthly_avg[,transcom_region_plot_monthly]
    new_data <- data.frame(MONTH =1:24, FLUX=c(subset_truth), KIND=rep("Truth",24))
    
    g = ggplot(combined_df, aes(x=MONTH, y= FLUX, fill=KIND)) +
      geom_boxplot(width=0.5) +   # outlier.shape = NA
      scale_fill_manual(values=c("red", "blue","green")) +
      geom_point(data=new_data, aes(x=MONTH, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
      labs(title=transcom_names[transcom_region_plot_monthly]) + theme(axis.text=element_text(size=12),
                                                                       axis.title=element_text(size=14,face="bold"),title=element_text(size=16),legend.text=element_text(size=14))
    
    
    
    print(g)
    
  }
  
  
  #-- Plot global sum 
  subset = ret$prior_tr_monthly_global
  prior_df = cbind(FLUX=as.vector(subset),KIND=rep("Prior",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
  
  subset = ret$post_tr_monthly_global
  post_df = cbind(FLUX=as.vector(subset),KIND=rep("Post",length(subset)),MONTH=rep(1:24,(dim(subset)[2])))
  
  combined_df  = as.data.frame(rbind(post_df,prior_df))
  
  combined_df$KIND = factor(combined_df$KIND)
  combined_df$FLUX = as.numeric(combined_df$FLUX)
  combined_df$MONTH = factor(combined_df$MONTH,levels=1:24)
  require(ggplot2)
  
  options(repr.plot.width=20, repr.plot.height=8)
  
  #-- This plot is monthly avg regional flux (in units of PgC/yr), plotted for 2 years
  
  subset_truth = ret$transcom_fluxes_real_monthly_avg_global
  new_data <- data.frame(MONTH =1:24, FLUX=c(subset_truth), KIND=rep("Truth",24))
  
  g = ggplot(combined_df, aes(x=MONTH, y= FLUX, fill=KIND)) +
    geom_boxplot(width=0.5) +   # outlier.shape = NA
    scale_fill_manual(values=c("red", "blue","green")) +
    geom_point(data=new_data, aes(x=MONTH, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
    labs(title="Global CO2 Flux") + theme(axis.text=element_text(size=12),
                                          axis.title=element_text(size=14,face="bold"),title=element_text(size=16),legend.text=element_text(size=14))
  
  
  
  print(g)
}

#-- Use this to control plot size
plot_timeseries_flux_bytranscom = function(data=ret)
{
  options(repr.plot.width=20, repr.plot.height=8)
 
  mat_table = data.frame(NUMBER=1:22,NAME=transcom_names)
  
  plt_df = ret$full_df
  plt_df$REGION = mat_table$NAME[plt_df$REGION]
                        
  #-- I believe this plot is monthly avg regional flux (in units of PgC/yr) for two years
  g = ggplot(plt_df, aes(x=REGION, y= FLUX, fill=KIND)) +
    geom_boxplot(width=0.5) +   # outlier.shape = NA
    scale_fill_manual(values=c("red", "blue","green"))
  
  
  new_data <- data.frame(REGION =transcom_names, FLUX=c(data$transcom_fluxes_real_annual_avg ), KIND=rep("Truth",22))
  
  g + geom_point(data=new_data, aes(x=REGION, y=FLUX, fill=KIND), color="black",bg="green", size=5, pch=21) +
    theme(axis.text.x = element_text(angle=70,size=15))
  
}

#-- FORM DATA FRAME OF PRIOR/POSTERIOR USING SAMPLES FROM NETCDF FILES, FOR PLOTTING
#-- this needs transparent parapply wrapper on the monthly grid2transcom() to speed up (e.g. line 2 below)
organize_data_for_plotting = function(prior_flux_netcdf="/Users/aschuh/test_output_prior/gridded_fluxes.nc4",
                                      posterior_flux_netcdf="/Users/aschuh/test_output_posterior/gridded_fluxes.nc4",max.cores=4)
{
  cores_used = min(max.cores,detectCores() - 2)
  cl <- makeCluster(cores_used, type = "FORK")
   print(paste("you have",detectCores(),"cores to work with, working with ",cores_used))
  d = load.ncdf(posterior_flux_netcdf)
  post_tr_monthly <- parApply(cl,d$flux.samples,c(3,4),FUN=function(x){grid2transcom(x,file_location=data_dir)}) %>% aperm(c(2,3,1)) 
  #post_tr_monthly = aaply(d$flux.samples,c(3,4),.fun=function(x){grid2transcom(x)})
  post_tr_monthly = post_tr_monthly[,,2:23]  #drop transcom=0
  post_tr_monthly_global = apply(post_tr_monthly,c(1,2),sum)
  post_tr_annual = apply(post_tr_monthly,c(2,3),sum)*30.5*3600*24*1e3*0.5 * 1e-15# ~ gC/yr
  post_df = cbind(FLUX=as.vector(t(post_tr_annual)),KIND=rep("Post",length(post_tr_annual)),REGION=rep(1:22,dim(post_tr_annual)[1]))
  
  e = load.ncdf(prior_flux_netcdf)
  prior_tr_monthly <- parApply(cl,e$flux.samples,c(3,4),FUN=function(x){grid2transcom(x,file_location=data_dir)}) %>% aperm(c(2,3,1))
  #prior_tr_monthly = aaply(e$flux.samples,c(3,4),.fun=function(x){grid2transcom(x)})
  prior_tr_monthly = prior_tr_monthly[,,2:23] 
  prior_tr_monthly_global = apply(prior_tr_monthly,c(1,2),sum)
  prior_tr_annual = apply(prior_tr_monthly,c(2,3),sum)*30.5*3600*24*1e3*0.5 * 1e-15 # ~ ~ PgC/yr adjustment to prior, 0.5 is 2 yr -> 1yr
  prior_df = cbind(FLUX=as.vector(t(prior_tr_annual)),KIND=rep("Prior",length(prior_tr_annual)),REGION=rep(1:22,dim(prior_tr_annual)[1]))
  
  full_df  = as.data.frame(rbind(post_df,prior_df))
  
  stopCluster(cl)
  
  full_df$KIND = factor(full_df$KIND)
  full_df$FLUX = as.numeric(full_df$FLUX)
  full_df$REGION = factor(full_df$REGION,levels=1:22)
  
  
  #-- This calculate the "true" fluxes from original prior fluxes
  transcom_fluxes_real = pull_true_transcom_flux(prior_flux_file=file.path(data_dir,"priors/prior_SiB4.nc"),state_true=state_vector_true)
  transcom_fluxes_real_annual_avg = transcom_fluxes_real$annual_2yr
  transcom_fluxes_real_monthly_avg = transcom_fluxes_real$monthly
  transcom_fluxes_real_monthly_avg_global = apply(transcom_fluxes_real$monthly,c(1),sum)
  
  output = list(prior_tr_monthly=prior_tr_monthly,prior_tr_monthly_global=prior_tr_monthly_global,prior_tr_annual=prior_tr_annual,
                post_tr_monthly=post_tr_monthly,post_tr_monthly_global=post_tr_monthly_global,post_tr_annual=post_tr_annual,
                full_df=full_df,transcom_fluxes_real_annual_avg=transcom_fluxes_real_annual_avg,
                transcom_fluxes_real_monthly_avg=transcom_fluxes_real_monthly_avg,
                transcom_fluxes_real_monthly_avg_global=transcom_fluxes_real_monthly_avg_global)
  return(output)
}

ar_covariance <- function(n, rho, variance = 1) {
  variance * toeplitz(rho ^ (0 : (n - 1)))
}

#-- Transcom region labels for regions 1:22
transcom_names = c("North American Boreal    ", "North American Temperate ",
                   "South American Tropical  ", "South American Temperate ",
                   "Northern Africa          ", "Southern Africa          ",
                   "Eurasia Boreal           ", "Eurasia Temperate        ",
                   "Tropical Asia            ", "Australia                ",
                   "Europe                   ", "North Pacific Temperate  ",
                   "West Pacific Tropical    ", "East Pacific Tropical    ",
                   "South Pacific Temperate  ", "Northern Ocean           ",
                   "North Atlantic Temperate ", "Atlantic Tropical        ",
                   "South Atlantic Temperate ", "Southern Ocean           ",
                   "Indian Tropical          ", "South Indian Temperate   ")

