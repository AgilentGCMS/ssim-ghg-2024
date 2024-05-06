#######################################################
#-- This code collects the GeosCHEM HEMCO flux output
#-- from Bertolacci Transcom case for Summer School 24'
#-- in order to create a "prior" which can be convoluted
#-- across inversion solution to create posterior flux est
#######################################################

setwd("")

#runs_folder = "runs_sib4_prior_4x5_mask"
runs_folder = "runs"

fls=list.files(paste("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/20240109_basis_functions/",runs_folder,sep=""),pattern="HEMCO_diagnostics",recursive=TRUE)
fls = fls[-grep("base",fls)]
#fls = fls[-grep("20140901_part001_split02",fls)]

fls_part001 = fls[grep("_part001",fls)]
fls_part002 = fls[grep("_part002",fls)]

mapp = read.csv(paste("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/20240109_basis_functions/",runs_folder,"/mapping.csv",sep=""))

require(ncdf4)
setwd(paste("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/20240109_basis_functions/",runs_folder,sep=""))

prior_part1 = array(0,dim=c(72,46,24))

for(i in 1:length(fls_part001)){
  print(paste("working on ",fls_part001[i]))
  con = nc_open(fls_part001[i])
  vrs = names(con$var)[grep("Emis",names(con$var))]
  for(j in 1:length(vrs))
  {
    print(j)
    id = mapp[ match(gsub("Emis_","",vrs[j]) , mapp[,3] ) , 1]
    print(id)
    if(substring(id,1,10) != "background"){
      datee = strsplit(id,"month")[[1]][2]
      yr = as.numeric(substr(datee,1,4) )
      mon = as.numeric(substr(datee,6,7) )
      val = ncvar_get(con,vrs[j],start=c(1,1,1,1),count=c(-1,-1,1,1))
      print(paste("index:",(yr-2014)*12+(mon-9) + 1))
      prior_part1[,,(yr-2014)*12+(mon-9) + 1] = prior_part1[,,(yr-2014)*12+(mon-9) + 1] + val
    }
  }
  nc_close(con)
}

prior_part2 = array(0,dim=c(72,46,24))

for(i in 1:length(fls_part002)){
  print(paste("working on ",fls_part002[i]))
  con = nc_open(fls_part002[i])
  vrs = names(con$var)[grep("Emis",names(con$var))]
  for(j in 1:length(vrs))
  {
    print(j)
    id = mapp[ match(gsub("Emis_","",vrs[j]) , mapp[,3] ) , 1]
    print(id)
    if(substring(id,1,10) != "background"){
      datee = strsplit(id,"month")[[1]][2]
      yr = as.numeric(substr(datee,1,4) )
      mon = as.numeric(substr(datee,6,7) )
      val = ncvar_get(con,vrs[j],start=c(1,1,1,1),count=c(-1,-1,1,1))
      print(paste("index:",(yr-2014)*12+(mon-9) + 1))
      prior_part2[,,(yr-2014)*12+(mon-9) + 1] = prior_part2[,,(yr-2014)*12+(mon-9) + 1] + val
    }
  }
  nc_close(con)
}

prior_part = prior_part1 + prior_part2 



t <- ncdim_def( "time", "months since 2014-09-01", 0:23, unlim=TRUE)

x <- ncdim_def( "longitude", units="degrees_east", vals=seq(-180,177.5,by=5))

y <- ncdim_def( "latitude", units="degrees_north", vals=c(-89,seq(-87,87,by=4),89))

newdim = list(x,y,t)  #,pft.dim)

varNEE <- ncvar_def(name="NEE",
                    units="kgCO2/m2/s",
                    dim=newdim,
                    missval=1e36,
                    longname="NEE Flux",
                    compression=5,
                    prec="float")

ncnew = nc_create(paste("/discover/nobackup/projects/csu/aschuh/temp/prior_SiB4_unit_pulse.nc",sep=""),
                  vars=list(varNEE),force_v4 = TRUE)

ncvar_put( ncnew, varNEE,prior_part,start=c(1,1,1),count=c(-1,-1,-1))

nc_close(ncnew)

#######################################################
#-- End of collection of the GeosCHEM HEMCO flux output
############################################################



###################################################################
#-- Use MIPv10 posterior flux results to create synthetic "truth"
#-- scalings and adjustments for toy example
#-- Using transport and prior from Sept 2014 to Aug 2016 and
#-- "truth" from Sept 2015 to Aug 2017
###################################################################

#-- Read SiB4 and OCO2 MIP posteriors and form differences and ratios for "truths" in toy example
require(ncdf4)
fil = nc_open("/discover/nobackup/projects/csu/aschuh/temp/prior_SiB4_unit_pulse.nc")
nee_prior = ncvar_get(fil,"NEE")
nc_close(fil)

#-- convert kgCO2/m2/sec to gC/m2/yr
nee_prior = nee_prior * 1e3 * 12/44 * 365 *24 * 3600

require(plyr)
ddir = "/discover/nobackup/projects/csu/aschuh/data/inversion_workshop_data/"

#-- this should output gC/yr/transcomregion
transcom_nee_prior = aaply(nee_prior,3,.fun=function(x){grid2transcom(x,model.grid.x=5,model.grid.y=4, file_location=ddir)})


#-- Pulling v10 MIP fluxes

v10_files = list.files("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/v10_flux_mip/gridded_fluxes/5x4/",full.names=TRUE)
v10_files_short = list.files("/discover/nobackup/projects/csu/aschuh/WOMBAT_stuff/v10_flux_mip/gridded_fluxes/5x4/",full.names=FALSE)

truth_array = array(NA,dim=c(24,23,length(v10_files),2))

for(i in 1:length(v10_files))
{
  print(paste("working on...",v10_files[i]))
  con = nc_open(v10_files[i])
  nee_post = ncvar_get(con,"land") + ncvar_get(con,"ocean")
  nee_post = nee_post[,,9:(9+23)]
  
  #-- this should output gC/yr/transcomregion
  transcom_nee_post = aaply(nee_post,3,.fun=function(x){grid2transcom(x,model.grid.x=5,model.grid.y=4, file_location=ddir)})
  
  #-- this should be ratio of total transcom region fluxes between post and prior
  truth_array[,,i,1] = transcom_nee_post/transcom_nee_prior-1
  
  #-- this should be difference of total transcom region fluxes between post and prior (gC per year)
  truth_array[,,length(v10_files),2] = transcom_nee_post - transcom_nee_prior
}

dimnames(truth_array)[[2]] =c("Not optimized","North American Boreal","North American Temperate","South American Tropical",
                              "South American Temperate ","Northern Africa" ,"Southern Africa",
                              "Eurasia Boreal" ,"Eurasia Temperate","Tropical Asia" ,
                              "Australia","Europe" ,"North Pacific Temperate",
                              "West Pacific Tropical" ,"East Pacific Tropical","South Pacific Temperate" ,
                              "Northern Ocean","North Atlantic Temperate" ,"Atlantic Tropical",
                              "South Atlantic Temperate" ,"Southern Ocean","Indian Tropical" ,
                              "South Indian Temperate")

dimnames(truth_array)[[1]] = c(paste("2014",sapply(9:12,pad,width=2,fill="0"),sep=""),
                               paste("2015",sapply(1:12,pad,width=2,fill="0"),sep=""),
                               paste("2016",sapply(1:8,pad,width=2,fill="0"),sep=""))

dimnames(truth_array)[[3]] = v10_files_short

dimnames(truth_array)[[4]] = c("Scaling","Difference")


#aa=vector()
#for(i in 1:90){#  a = apply(truth_array[,2:23,i,1] * transcom_nee_prior[,2:23],1,sum,na.rm=TRUE)
#  b = apply(tm(truth_array[,2:23,i,1],-1,1) * transcom_nee_prior[,2:23],1,sum,na.rm=TRUE)
#  aa[i] = round(sum((b-a)/a),3)
#}

#save(truth_array,file="/discover/nobackup/projects/csu/aschuh/temp/truth_array_unit_pulse.rda")

######################################################
#-- END of scalings and adjustments for toy example
######################################################