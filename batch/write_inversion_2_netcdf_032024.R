write_inversion_2_netcdfs=function(inv_object=ret2,prior_mean_ncdf=file.path(data_dir,"priors/prior_SiB4.nc"),
                          #sensitivity.matrix="/discover/nobackup/projects/csu/aschuh/data/inversion_workshop_data/rdas/trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda",
                          sample_number=40,output_dir=output_dir)
{
#-- set.seed ?  
  
#-- Could be mismatch between inv_object's sens matrix and what is loaded here
#-- need to fix that if we pull in sens matrix

#-- required libs
  
require(ncdf4)
require(mvtnorm)
require(plyr)  
require(dplyr)
  
if(dir.exists(output_dir)){
  print(paste("overwriting inversion results in",output_dir))
}else{
  print(paste("creating",output_dir))
  dir.create(output_dir,recursive=TRUE)
}

outfile_name_state_vector = file.path(output_dir,"state_vector.nc4")
outfile_name_gridded_fluxes = file.path(output_dir,"gridded_fluxes.nc4")
outfile_name_cosamples = file.path(output_dir,"cosamples.nc4")
#load(sensitivity.matrix)
#obs_len = dim(jacob)[[1]]
  
  
####################################
#--  create dims for all netcdf vars
####################################

t <- ncdim_def( "time", "months since 2014-09-01", 0:23, unlim=TRUE)

state_dim <- ncdim_def( "scaling_factor_number", "unitless", 1:528, unlim=TRUE)

sample_dim <- ncdim_def( "sample_dim", "", 1:sample_number, unlim=TRUE)

x_dim <- ncdim_def( "longitude", units="degrees_east", vals=seq(-179.5,179.5,by=1))

y_dim <- ncdim_def( "latitude", units="degrees_north", vals=seq(-89.5,89.5,1))

cosample_dim <- ncdim_def( "cosamples", units="", vals=1:length(inv_object$inputs$obs_errors))

#######################################
#--  create netcdf vars for statevector
#######################################

varSTATE_MEAN = ncvar_def(name="scaling_factor_mean",
                         units="",
                         dim=state_dim,
                         prec="float")

varSTATE_SAMPLE = ncvar_def(name="scaling_factor_sample",
                          units="",
                          dim=list(state_dim,sample_dim),
                          prec="float")

varSTATE_NAME = ncvar_def(name="scaling_factor_name",
          units="",
          dim=state_dim,
          prec="char")

varSTATE_COV = ncvar_def(name="scaling_factor_covariance",
                        units="",
                        dim=list(state_dim,state_dim),
                        prec="float")

#######################################
#--  create netcdf vars for cosamples
#######################################

varOBS_TYPE = ncvar_def(name="observation_type",
                         units="mole fraction CO2",
                         dim=cosample_dim,
                         prec="char")

varOBS_ID = ncvar_def(name="observation_id",
                       units="mole fraction CO2",
                       dim=cosample_dim,
                       prec="char")

varOBS_MEAN = ncvar_def(name="mole_fraction_mean",
                       units="mole fraction CO2",
                       dim=cosample_dim,
                       prec="float")

varOBS_SD = ncvar_def(name="mole_fraction_sd",
                       units="mole fraction CO2",
                       dim=cosample_dim,
                       prec="float")

varOBS_ASSIM_FLAG = ncvar_def(name="assim flag",
                      units="mole fraction CO2",
                      longname="assim,0=NO,1=YES",
                      dim=cosample_dim,
                      prec="integer")

varOBS_INPUT = ncvar_def(name="original input obs",
                         units="mole fraction CO2",
                         dim=cosample_dim,
                         prec="float")
  
##########################################
#--  create netcdf vars for gridded fluxes
##########################################

varNEE_MEAN <- ncvar_def(name="flux_mean",
                    units="kgCO2/m2/s",
                    dim=list(x_dim,y_dim,t),
                    missval=1e36,
                    longname="NEE Flux",
                    compression=5,
                    prec="float")

varNEE_SD <- ncvar_def(name="flux_sd_marginal",
                         units="kgCO2/m2/s",
                         dim=list(x_dim,y_dim,t),
                         missval=1e36,
                         longname="NEE Flux SD",
                         compression=5,
                         prec="float")

varNEE_samples <- ncvar_def(name="flux_samples",
                       units="kgCO2/m2/s",
                       dim=list(x_dim,y_dim,t,sample_dim),
                       missval=1e36,
                       longname="NEE Flux Sample",
                       compression=5,
                       prec="float")

##########################################
#--  pull samples for output
##########################################

print("pulling multivariate normal sample....")

samps = rmvnorm(sample_number,mean = inv_object$x_hat,sigma=inv_object$P)

samps_array = aaply(samps,1,.fun=function(x){matrix(x,nrow=24,byrow=FALSE)})
  
##########################################
#--  create state_vector netcdf
##########################################
print(paste("writing state vector output..to ",outfile_name_state_vector))

if(file.exists(outfile_name_state_vector)){file.remove(outfile_name_state_vector)}

ncnew_statevector = nc_create(outfile_name_state_vector,
                  vars=list(varSTATE_MEAN,varSTATE_SAMPLE,varSTATE_NAME,varSTATE_COV),force_v4 = TRUE)

ncvar_put(ncnew_statevector,varSTATE_MEAN,vals=inv_object$x_hat)

ncvar_put(ncnew_statevector,varSTATE_SAMPLE,vals=t(samps))

ncvar_put(ncnew_statevector,varSTATE_NAME,vals=dimnames(jacob)[[2]])

ncvar_put(ncnew_statevector,varSTATE_COV,vals=inv_object$P)
          
nc_close(ncnew_statevector)

##########################################
#--  create 1x1 gridded_fluxes for netcdf
##########################################
print("creating gridded fluxes....")

con = nc_open(prior_mean_ncdf)
NEE = ncvar_get(con,"NEE")
nc_close(con)
NEE_1x1 = aaply(NEE,3,.fun=function(x){expand_5x4_2_1x1(x)}) %>%
           aperm(c(2,3,1))

tr_dir = file.path(data_dir,"/transcom/",sep="")

x_hat_matrix = matrix(inv_object$x_hat,nrow=24,byrow=FALSE)

gridded_1x1_state = aaply(x_hat_matrix,1,.fun=function(x){
             transcom2grid(x,model.grid.x = 1,model.grid.y = 1,
              file_location = data_dir)}) %>% aperm(c(2,3,1))

gridded_1x1_mean_flux = NEE_1x1*gridded_1x1_state

gridded_1x1_state_samples = aaply(samps_array,c(1,2),.fun=function(x){transcom2grid(x,model.grid.x = 1,
                                      model.grid.y = 1,file_location = data_dir)}) %>%
                             aperm(c(3,4,2,1))

gridded_1x1_flux_samples = array(NA,dim=c(dim(gridded_1x1_mean_flux),sample_number))

    
#-- need to look for cores to use and foreach/parallel/parapply
#for(i in 1:sample_number){
#  tmp_state = gridded_1x1_state_samples[,,i]
#  gridded_1x1_flux_samples[,,,i] = aaply(gridded_1x1_mean_flux,3,.fun=function(x){x*tmp_state}) %>%
#    aperm(c(2,3,1))
#}

gridded_1x1_flux_samples = aaply(gridded_1x1_state_samples,c(4),.fun=function(x){ x*NEE_1x1 }) %>% aperm(c(2,3,4,1))
                 

gridded_1x1_flux_sd = apply(gridded_1x1_flux_samples,c(1,2,3),sd)


##########################################
#--  write gridded_fluxes to netcdf
##########################################

print(paste("writing gridded fluxes...to ",outfile_name_gridded_fluxes))

if(file.exists(outfile_name_gridded_fluxes)){file.remove(outfile_name_gridded_fluxes)}

ncnew_gridded_fluxes = nc_create(outfile_name_gridded_fluxes,
                     vars=list(varNEE_MEAN,varNEE_SD,varNEE_samples),force_v4 = TRUE)

ncvar_put(ncnew_gridded_fluxes,varNEE_MEAN,vals=gridded_1x1_mean_flux)

ncvar_put(ncnew_gridded_fluxes,varNEE_samples,vals=gridded_1x1_flux_samples)

ncvar_put(ncnew_gridded_fluxes,varNEE_SD,vals=gridded_1x1_flux_sd)

nc_close(ncnew_gridded_fluxes)

##########################################
#--  collect cosamples output
##########################################

print("creating cosamples output....")

#obs_mean = jacob %*% inv_object$x_hat
obs_mean = inv_object$outputs$modeled_obs

#-- What are we talking about here?  Obs errors in inversion or posterior sd around model estimate of obs?
obs_sd = inv_object$inputs$obs_errors

obs_id = inv_object$inputs$obs_id

obs_orig = inv_object$inputs$obs

obs_type = rep("SAT_XCO2",dim(jacob)[1])
obs_type[grep("obspack",obs_id)] = "INSITU"
nch = sapply(obs_id,nchar)
obs_type[nch==14] = "TCCON_XCO2"

obs_assim = inv_object$inputs$assim_T_F

##########################################
#--  create cosamples netcdf
##########################################

print(paste("writing cosamples output...to ",outfile_name_cosamples))

if(file.exists(outfile_name_cosamples)){file.remove(outfile_name_cosamples)}

ncnew_cosamples = nc_create(outfile_name_cosamples,
                                 vars=list(varOBS_TYPE,varOBS_ID,varOBS_MEAN,varOBS_SD,varOBS_ASSIM_FLAG,varOBS_INPUT),force_v4 = TRUE)
 
ncvar_put(ncnew_cosamples,varOBS_TYPE,vals=obs_type)

ncvar_put(ncnew_cosamples,varOBS_ID,vals=obs_id)

ncvar_put(ncnew_cosamples,varOBS_MEAN,vals=obs_mean)

ncvar_put(ncnew_cosamples,varOBS_INPUT,vals=obs_orig)

ncvar_put(ncnew_cosamples,varOBS_SD,vals=obs_sd)

ncvar_put(ncnew_cosamples,varOBS_ASSIM_FLAG,vals=obs_assim)

nc_close(ncnew_cosamples)
}


#write_inversion_2_netcdfs(inv_object=ret2)