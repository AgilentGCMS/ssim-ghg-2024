#######################################################
#-- ***Parent Directory and code for ALL inversions***
#######################################################

#dir = "/discover/nobackup/projects/csu/aschuh/data/inversion_workshop_data/"
#dir = "/discover/nobackup/projects/csu/aschuh/SSIM-GHG/base"
#dir = "/Users/aschuh/SSIM-GHG/base"
dir = "../"

###############################################
#-- Load Code
##############################################
source(file.path(dir,"R_code/util_code_031024.R"))
source(file.path(dir,"R_code/inversion_031024.R"))
source(file.path(dir,"R_code/write_inversion_2_netcdf_031024.R"))
###############################################
#--  Load sensitivity matrices 
###############################################

load(file.path(dir,"data/","trunc_full_jacob_030624_with_dimnames_sib4_4x5_mask.rda"))
load(file.path(dir,"data/","jacob_bgd_021624.rda"))


###################################################################
#-- END END END ***Parent Directory and code for ALL inversions***
###################################################################

##################################################################
#- Inversion #1   *************************
##################################################################

#################################
#- Target truth in state space
#################################

##################################################################
#-- This array holds ratios of OCO2v10MIP fluxes and SiB4 fluxes
#-- as examples of "scalings" to be recovered. It also holds corresponding
#-- differences if the inversion attempts to directly solve for flux
##################################################################

#load("/projects/sandbox/inversion_workshop_scripts/truth_array.rda")
load(file.path(dir,"data/truth_array.rda"))

xx = truth_array[,-1,1,1]

state_vector_true= tm(as.vector(truth_array[,-1,1,1]),-1,1)

#########################################################
# Generate a prior flux covariance matrix P_0
# Long term, a catalog of predefined choices is best here I think
#########################################################
land_prior_sd = 0.3
ocean_prior_sd = 0.3
require(Matrix)
require(MixMatrix)

#-- induce temporal correlations
temporal_cov = list(ARgenerate(n=24, rho=0.5))
sigma = bdiag(rep(temporal_cov,11*2))

#-- scale by variance for land/ocean
var_scaling_diagonal = diag(c(rep(land_prior_sd,24*11),rep(ocean_prior_sd,24*11)))
sigma = var_scaling_diagonal %*% sigma %*% t(var_scaling_diagonal)

##########################################################
#-- sd for Gaussian i.i.d. errors, jacob is sens matrix
##########################################################
R_diagonal_in = rep(1,(dim(jacob)[1]))

##########################################
#-- Generate obs, 'y',  set.seed() ????
##########################################

y_in = generate_observations(H=jacob,H_bgd=jacob_bgd,
                             state_vector=state_vector_true,err_obs=R_diagonal_in)

####################################################################################
#-- WHICH obs do you want to use in the inversion? 
#-- examples of selecting on stations, type of data, lat/lon box,etc
####################################################################################

load(file.path(dir,"data/obs_catalog_030624.rda")) # obs_catalog object

subset_indicator_obs=rep(FALSE,dim(jacob)[1])

############################
#-- SAMPLE BY TYPE EXAMPLE
############################
#subset_indicator_obs[obs_catalog$TYPE == "TCCON"] = TRUE

############################
#-- SAMPLE BY NOAA STATION EXAMPLE
############################
#subset_indicator_obs[grep("spo",obs_catalog$ID)] = TRUE
#subset_indicator_obs[grep("lef",obs_catalog$ID)] = TRUE

############################
#-- SAMPLE BY TIME EXAMPLE
############################
#subset_indicator_obs[obs_catalog$TIME < 8000000] = TRUE

############################
#-- SAMPLE BY LON & LAT EXAMPLE
############################
#subset_indicator_obs[obs_catalog$LON < -10 & obs_catalog$LAT > 10] = TRUE

subset_indicator_obs=c(rep(TRUE,156383),rep(FALSE,1000000))

############################
#-- Run the actual inversion
############################

ret2 = invert_clean(H=jacob,R_diagonal=R_diagonal_in,P_0=sigma,y=y_in,H_bgd=jacob_bgd,
                    subset_indicator_obs=c(rep(TRUE,156383),rep(FALSE,1000000)))


#############################################
#-- Output data from inversion to netcdf
#############################################

#write_inversion_2_netcdfs(inv_object=ret2,output_dir="/discover/nobackup/projects/csu/aschuh/temp/example_inversion")
write_inversion_2_netcdfs(inv_object=ret2,
                          prior_mean_ncdf=file.path(dir,"data/prior_SiB4.nc"),
                          output_dir="/Users/aschuh/test_output")

sessionInfo()