generate_transcom_flux_ensemble_from_inversion=function(inv_object=ret2,
                                                        prior_mean_ncdf=file.path(data_dir,"priors/prior_SiB4.nc"),
                                                        samples=500)
{
  #-- set.seed ?  
  
  require(ncdf4)
  require(mvtnorm)
  require(plyr)  
  require(dplyr)
  
  con = nc_open(prior_mean_ncdf)
  NEE_5x4 = ncvar_get(con,"NEE")
  nc_close(con)
  
  #NEE_1x1 = aaply(NEE,3,.fun=function(x){expand_5x4_2_1x1(x)}) %>%
  #  aperm(c(2,3,1))
  
  x_hat_matrix = matrix(inv_object$posterior$x_hat,nrow=24,byrow=FALSE)
  
  tr_dir = file.path(data_dir,"/transcom/",sep="")
  
  transcom_NEE_5x4 = aaply(NEE_5x4,c(3),.fun=function(x){grid2transcom(x,model.grid.x = 5,model.grid.y = 4)})
  transcom_NEE_5x4 = transcom_NEE_5x4[,-1]
  
  x_hat_matrix_v = as.vector(x_hat_matrix)
  transcom_NEE_5x4_v = as.vector(transcom_NEE_5x4)
  
  transcom_NEE_5x4_reals = rmvnorm(n=samples,mean=transcom_NEE_5x4_v*x_hat_matrix_v,sigma = (diag(transcom_NEE_5x4_v) %*% ret2$posterior$Sx_post %*% diag(transcom_NEE_5x4_v))) %>% aperm(c(2,1))
  transcom_NEE_5x4_reals = aaply(transcom_NEE_5x4_reals,c(2),.fun=function(x){array(x,dim=c(24,22))}) %>% aperm(c(2,1,3))
  
  post_tr_monthly = transcom_NEE_5x4_reals     *12/44    *30.5*3600*24*1e3 * 1e-15
  post_tr_monthly_global = apply(post_tr_monthly,c(1,2),sum)
  #post_tr_annual = apply(post_tr_monthly,c(2,3),sum)  *12/44   *30.5*3600*24*1e3*0.5 * 1e-15# ~ gC/yr
  post_tr_annual = apply(post_tr_monthly,c(2,3),sum) * 0.5 # ~ PgC/yr
  post_df = cbind(FLUX=as.vector(t(post_tr_annual)),KIND=rep("Post",length(post_tr_annual)),REGION=rep(1:22,dim(post_tr_annual)[1]))
  
  transcom_NEE_5x4_reals = rmvnorm(n=samples,mean=rep(0,length(x_hat_matrix_v)),sigma = (diag(transcom_NEE_5x4_v) %*% ret2$prior$Sx %*% diag(transcom_NEE_5x4_v))) %>% aperm(c(2,1))
  transcom_NEE_5x4_reals = aaply(transcom_NEE_5x4_reals,c(2),.fun=function(x){array(x,dim=c(24,22))}) %>% aperm(c(2,1,3))
  
  prior_tr_monthly = transcom_NEE_5x4_reals       *  12/44 *30.5*3600*24*1e3* 1e-15
  prior_tr_monthly_global = apply(prior_tr_monthly,c(1,2),sum)
  #prior_tr_annual = apply(prior_tr_monthly,c(2,3),sum)  * 12/44   *30.5*3600*24*1e3*0.5 * 1e-15 # ~ ~ PgC/yr adjustment to prior, 0.5 is 2 yr -> 1yr
  prior_tr_annual = apply(prior_tr_monthly,c(2,3),sum)  * 0.5  # ~ ~ PgC/yr adjustment to prior, 0.5 is 2 yr -> 1yr
  prior_df = cbind(FLUX=as.vector(t(prior_tr_annual)),KIND=rep("Prior",length(prior_tr_annual)),REGION=rep(1:22,dim(prior_tr_annual)[1]))
  
  full_df  = as.data.frame(rbind(post_df,prior_df))
  
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
  
