#################################################################################################
#-- PART I:  Generate Jacobians from offline sensitivities (from GEOS-Chem pulses)
#################################################################################################

#-- Each directory effectively contains one column of the Obs by State Jacobian
dir_list = list.dirs("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/20240109_basis_functions/sensitivities_unit_flux_4x5_mask/")
base_inds = grep("base_",dir_list)
bgds = dir_list[base_inds]
dir_list = dir_list[-c(1,base_inds)] # first dir appears to be root directory
#transcom_ocn_drops = sapply(1:11,FUN=function(x){grep(paste("ocean_regionRegion",pad(x,width=2,fill="0"),sep=""),dir_list)})
#transcom_land_drops = sapply(12:22,FUN=function(x){grep(paste("nee_regionRegion",pad(x,width=2,fill="0"),sep=""),dir_list)})
#dir_list = dir_list[-c(transcom_ocn_drops,transcom_land_drops)]
#dir_nee_list = dir_list[grep("nee_",dir_list)]
#dir_ocn_list = dir_list[grep("ocean_",dir_list)]

#-- Generate the basic Jacobian for inversion from sensitivities
for(i in 1:length(dir_list))
{
  print(paste("reading",dir_list[i]))
  oco2_data = load.ncdf(paste(dir_list[i],"/oco2-sensitivity.nc4",sep=""))
  tccon_data = load.ncdf(paste(dir_list[i],"/tccon-sensitivity.nc4",sep=""))
  obspack_data = load.ncdf(paste(dir_list[i],"/obspack-assim-1-sensitivity.nc4",sep=""))
  if(i==1){
    # Double check there are no duplicates
    stopifnot(
      length(obspack_data$observation.id) == length(unique(as.character(obspack_data$observation.id))),
      length(tccon_data$observation.id) == length(unique(as.character(tccon_data$observation.id))),
      length(oco2_data$observation.id) == length(unique(as.character(oco2_data$observation.id))),
      length(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)) == 
        length(unique(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)))
    )
    #-- build observation "catalog" for filtering
    all_observation_ids <- c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)
    all_observations_ids_save <- all_observation_ids
    
    all_time <- c(oco2_data$model.time,tccon_data$model.time,obspack_data$model.time)
    all_time_save <- all_time
    
    all_lat <- c(oco2_data$model.latitude,tccon_data$model.latitude,obspack_data$model.latitude)
    all_lat_save <- all_lat
    
    all_lon <- c(oco2_data$model.longitude,tccon_data$model.longitude,obspack_data$model.longitude)
    all_lon_save <- all_lon
    
    obs_type = c(rep("OCO2",length(oco2_data$model.longitude)),
                 rep("TCCON",length(tccon_data$model.longitude)),
                 rep("IS",length(obspack_data$model.longitude)))
    
    jacob <- array(0, dim = c(length(all_observation_ids), 22 * 24 * 2))
  }
  jacob[
    match(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id), 
          all_observation_ids),
    i
  ] = c(oco2_data$sensitivity,tccon_data$sensitivity,obspack_data$sensitivity)
}

#-- Give meaningful dimension names to Jacobian
nms = sapply(dir_list,FUN=function(x){strsplit(x,"//")[[1]][2]})
dimnames(jacob)[[2]] = as.character(nms)
dimnames(jacob)[[1]] = all_observations_ids_save

obs_catalog = data.frame(TYPE=obs_type,ID=all_observations_ids_save,LON=all_lon_save,LAT=all_lat_save,TIME=all_time_save)

#################################################################################################
#save(jacob,file="/discover/nobackup/projects/csu/aschuh/temp/full_jacob_030624_with_dimnames_sib4_4x5_mask.rda")
#jacob=jacob[,c(1:264,793:1056)];save(jacob,file="/discover/nobackup/projects/csu/aschuh/temp/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda")
#save(jacob,file="/discover/nobackup/projects/csu/aschuh/temp/full_jacob_032624_with_dimnames_unit_pulse_4x5_mask_hour_timestamping.rda")
#jacob=jacob[,c(1:264,793:1056)];save(jacob,file="/discover/nobackup/projects/csu/aschuh/temp/trunc_full_jacob_032624_with_dimnames_unit_pulse_4x5_mask_hour_timestamping.rda")
#save(obs_catalog,file="/discover/nobackup/projects/csu/aschuh/temp/obs_catalog_040324_unit_pulse_hour_timestamp.rda")
#################################################################################################

#-- Generate the basic Jacobian for fixed contributions, e.g. land prior, ocn prior, fire, and fossil, from sensitivities
#-- BE AWARE, land prior and ocn prior won't perfectly match H %*% 1 since the serial runs grab hourly output and the Jacobians
#-- grab hourly for first 4 months and daily avg afterwards

#-- Each directory effectively contains one column of the Obs by State Jacobian
dir_list = list.dirs("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/20240109_basis_functions/sensitivities_sib4_prior/")
base_inds = grep("base",dir_list)
bgds = dir_list[base_inds]
bgds = bgds[c(2,3,4,5,1)]

for(i in 1:length(bgds))
{
  print(paste("reading",bgds[i]))
  oco2_data = load.ncdf(paste(bgds[i],"/oco2-sensitivity.nc4",sep=""))
  tccon_data = load.ncdf(paste(bgds[i],"/tccon-sensitivity.nc4",sep=""))
  obspack_data = load.ncdf(paste(bgds[i],"/obspack-assim-1-sensitivity.nc4",sep=""))
  if(i==1){
    # Double check there are no duplicates
    stopifnot(
      length(obspack_data$observation.id) == length(unique(as.character(obspack_data$observation.id))),
      length(tccon_data$observation.id) == length(unique(as.character(tccon_data$observation.id))),
      length(oco2_data$observation.id) == length(unique(as.character(oco2_data$observation.id))),
      length(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)) == 
        length(unique(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)))
    )
    all_observation_ids <- c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)
    all_observations_ids_save <- all_observation_ids
    jacob_bgd <- array(0, dim = c(length(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id)), 5))
  }
  #-- Match to above full Jacobian obs_id_list
  jacob_bgd[
    match(c(oco2_data$observation.id,tccon_data$observation.id,obspack_data$observation.id), 
          all_observations_ids_save),
    i
  ] = c(oco2_data$sensitivity,tccon_data$sensitivity,obspack_data$sensitivity)
}

#-- Give meaningful dimension names to Jacobian
nms = sapply(bgds,FUN=function(x){strsplit(x,"//")[[1]][2]})
dimnames(jacob_bgd)[[2]] = as.character(nms)
dimnames(jacob_bgd)[[1]] = all_observations_ids_save

#################################################################################################
#save(jacob_bgd,file="/discover/nobackup/projects/csu/aschuh/temp/jacob_bgd_050924.rda")
#################################################################################################

# Sep 2019: operation_mode has been superceded by data_type
# This is not the usual operation mode, it conflates surface type
# with operation mode.  DFB has coded it into the last digit of the
# sounding_id.

# 1=land nadir, 2=land glint, 3=land target, 4=land transition,
# 5=ocean nadir, 6=ocean glint, 7=ocean target, and 8=ocean
# transition 

#- MODIFY OBS_CATALOG TO ADD A FEW THINGS
load("/projects/SSIM-GHG/data/obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors.rda")
obs_catalog$DATE = as.POSIXct(obs_catalog$TIME*60,origin = "2000-01-01",tz = "GMT")
obs_catalog$YEAR = year(obs_catalog$DATE)
obs_catalog$MONTH = month(obs_catalog$DATE)
obs_catalog$SITE = NA
#ind_IS = obs_catalog$TYPE == "IS"
#obs_catalog[ind_IS,"SITE"] = sapply(obs_catalog[ind_IS,"ID"],
#                                    FUN=function(x){strsplit(strsplit(x,"~")[[1]][2],"-")[[1]][1] })
save(obs_catalog,file="/projects/SSIM-GHG/data/obs/obs_catalog_042424_unit_pulse_hour_timestamp_witherrors_withdates.rda")