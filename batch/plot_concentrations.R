plot_concentrations = function(inversion=NULL,prior_netcdf="~/temp/output/test_output_prior/cosamples.nc4",
                               posterior_netcdf="~/temp/output/test_output_posterior/cosamples.nc4",
                               site_strings=c("brw","mlo","smo","co2_spo_surface-flask","lef","wkt","wbi","nwr","hun"),
                               ggplotly=FALSE,add_prior_nee=FALSE,add_fossil=FALSE,add_fire=FALSE)
{
if(!is.null(inversion))
{
  ret2 = inversion
  
  if(!add_prior_nee){pr = - ret2$prior$outputs$modeled_obs}else{pr = rep(0,length(ret2$prior$outputs$modeled_obs))}
  if(add_fossil){pr = pr+ fossil_fixed}  else{}
  if(add_fire){pr = pr+ fire_fixed}  else{}
  df = data.frame(ID=rep(obs_catalog$ID,3),
                  DATE=rep(obs_catalog$DATE,3),
                  VALUE=c(ret2$posterior$inputs$obs + pr,
                          ret2$prior$outputs$modeled_obs + pr,
                          ret2$posterior$outputs$modeled_obs + pr),
                  TYPE=c(rep("OBS",dim(obs_catalog)[1]),
                         rep("PRIOR",dim(obs_catalog)[1]),
                         rep("POSTERIOR",dim(obs_catalog)[1])))
  
  df$TYPE =factor(df$TYPE)
  df$VALUE = as.numeric(df$VALUE)
  
}else{
  cosamples_prior = load.ncdf(prior_netcdf)
  cosamples_post = load.ncdf(posterior_netcdf)
  
  #df = data.frame(ID=rep(obs_catalog$ID,3),
  #                DATE=rep(obs_catalog$DATE,3),
  #                VALUE=c(ret2$posterior$inputs$obs,
  #                        ret2$prior$outputs$modeled_obs,
  #                        ret2$posterior$outputs$modeled_obs),
  #                TYPE=c(rep("OBS",dim(obs_catalog)[1]),
  #                       rep("PRIOR",dim(obs_catalog)[1]),
  #                       rep("POSTERIOR",dim(obs_catalog)[1])))
}
  
for(i in 1:length(site_strings))
{
  
 ind = 1:length(df$ID) %in% grep(site_strings[i],df$ID)
  
 df2 = df[ind,]
  
 options(repr.plot.width=20, repr.plot.height=8)
 
  g =  ggplot(df2, aes(x = DATE, y = VALUE,color=TYPE)) + 
    scale_color_manual(values = c("blue","red","black")) +
    geom_point(size = 0.05) + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 12), se = FALSE) +
    labs(title=paste(site_strings[i],"Time series"),x ="Date", y = "CO2") + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"),
          title=element_text(size=16),legend.text=element_text(size=14))
  
  if(ggplotly){ggplotly(g)}else{print(g)}
 }
}
  
  
  