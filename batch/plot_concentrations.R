plot_concentrations = function(inversion=NULL,prior_netcdf="~/temp/output/test_output_prior/cosamples.nc4",
                               posterior_netcdf="~/temp/output/test_output_posterior/cosamples.nc4",
                               site_strings=c("brw","mlo","smo","spo","lef","wkt","wbi","nwr","hun"),
                               ggplotly=FALSE)
{
library(plotly)
if(!is.null(inversion))
{
  ret2 = inversion
  
  df = data.frame(ID=rep(obs_catalog$ID,3),
                  DATE=rep(obs_catalog$DATE,3),
                  VALUE=c(ret2$posterior$inputs$obs,
                          ret2$prior$outputs$modeled_obs,
                          ret2$posterior$outputs$modeled_obs),
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
  
  g =  ggplot(df2, aes(x = DATE, y = VALUE,color=TYPE)) + 
    scale_color_manual(values = c("blue","red","black")) +
    geom_point(size = 0.05) + 
    geom_smooth(method = "lm", formula = y ~ poly(x, 12), se = FALSE) +
    labs(title=paste(site_strings[i],"Time series"),x ="Date", y = "CO2")
  
  if(ggplotly){ggplotly(g)}else{print(g)}
 }
}
  
  
  
  
  
  
  
cosamples_prior = load.ncdf("~/temp/output/test_output_prior/cosamples.nc4")
cosamples_post = load.ncdf("~/temp/output/test_output_posterior/cosamples.nc4")

dat = cbind(cosamples_prior$observation.id,cosamples_prior$mole.fraction.mean,
            cosamples_post$mole.fraction.mean,cosamples_post[["assim flag"]])

dat = 

df <- economics %>% dplyr::select(date, psavert, uempmed) %>%
  gather(key = "variable", value = "value", -date)
head(df, 3)
# Multiple line plot
g = ggplot(df, aes(x = date, y = value)) + 
  geom_line(aes(color = variable), size = 1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) +
  theme_minimal()

ggplotly(g)