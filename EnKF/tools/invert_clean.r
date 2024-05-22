invert_clean = function(H,R_diagonal,P_0,y,H_bgd,subset_indicator_obs=NULL, DOF=FALSE)
{
  # H=jacob;R_diagonal=R_diagonal_in;P_0=sigma;y=y_in;H_bgd=jacob_bgd

  #-- capture full sensitivity matrix info before subsetting
  if(is.null(subset_indicator_obs)){subset_indicator_obs=rep(TRUE,length(R_diagonal))}

  obs_id = dimnames(H)[[1]]
    obs_len_init = dim(H)[1]
      assim_T_F=rep(0,dim(H)[1])
        assim_T_F[subset_indicator_obs] = 1
	  state_names = dimnames(H)[[2]]

  inputs_list = list(obs_errors=R_diagonal,assim_T_F=assim_T_F,obs_id=obs_id,obs=y)


  #-- Initial subsetting of sensitivity matrix and observations
    #-- Note sensitivity matrix H is now *subset*
        H_non_assimilated=H[!subset_indicator_obs,]
	    H_bgd_non_assimilated=H_bgd[!subset_indicator_obs,]
	        H=H[subset_indicator_obs,]
		    H_bgd=H_bgd[subset_indicator_obs,]
		        R_diagonal = R_diagonal[subset_indicator_obs]
			    y = y[subset_indicator_obs]
			      #}

  #-- CORE INVERSION CODE BELOW

  R_diagonal_inv = R_diagonal^-1
    R2_diagonal_inv = R_diagonal^-2

  HR = H * R_diagonal_inv

  print("...cross product")

  tH_H = crossprod(HR)

  rm(HR)
    gc()

  computation_1 = diag(diag(P_0)^-1) + tH_H

  rm(tH_H)

  print(".. .deriving posterior covariance matrix of state, P")
    P = solve(computation_1)

  #---  Digression to calc S, the influence/sensivity matrix, and associated DOF
    #---  need: diag(R2_diagonal_inv[subset_indicator_obs])%*%(H[subset_indicator_obs,])%*%P%*%t(H[subset_indicator_obs,])
      #---  calculated as tr(X %*% Y) as sum across rows of X * t(Y)

if(DOF){
  S1 = R2_diagonal_inv * (H)
    S2 = t(P%*%t(H))
      S3 = S1*S2
        trS = apply(S3,1,sum)
	  sum_trS = sum(trS)
	    sum_trB = dim(P)[1] - sum(trS)
	      print(paste("DFS:",sum_trS))
	        print(paste("DFB:",sum_trB))
		  #----
		    #----
		      }else{sum_trS=NULL;sum_trB=NULL}

  apriori_obs = H %*% rep(1+0,dim(H)[2])
                          #generate_observations(H=H,H_bgd=H_bgd,
			  #state_vector=rep(0,(dim(H)[2])),err_obs=NULL)

  if(sum(!assim_T_F) == 1){
        apriori_obs_non_assimilated = sum(H_non_assimilated * rep(1+0,length(H_non_assimilated))) #+ sum(H_bgd_non_assimilated[2:3])
	    }else{
	          apriori_obs_non_assimilated = H_non_assimilated %*% rep(1+0,dim(H_non_assimilated)[2]) #+ sum(H_bgd_non_assimilated[2:3])

#generate_observations(H=H_non_assimilated,H_bgd=H_bgd_non_assimilated,
#state_vector=rep(0,(dim(H_non_assimilated)[2])),
                                             #err_obs=NULL)
					         }


  computation_3 = matrix(y-(apriori_obs),ncol=1)

  HR2 = H * R2_diagonal_inv

  HR2 = t(HR2)

  computation_4 = HR2  %*% computation_3 #+ xtra

  rm(HR2)

  print("...deriving posterior mean state, X_hat")
    x_hat = P %*% computation_4

  rm(computation_4)
    rm(computation_3)

#-- CORE INVERSION CODE ABOVE

#-- Construct the modeled observations (obs assimilated + obs nonassimilated)

  assim_obs =  (H) %*% x_hat+ (H %*% rep(1,dim(H)[2]))
    non_assim_obs =  (H_non_assimilated) %*% x_hat + (H_non_assimilated %*% rep(1,dim(H_non_assimilated)[2]))
      modeled_obs = rep(NA,obs_len_init)
        modeled_obs[subset_indicator_obs] = assim_obs
	  modeled_obs[!subset_indicator_obs] = non_assim_obs

  #-- Construct the a priori observations
    apriori_obs_out = rep(NA,obs_len_init)
      apriori_obs_out[subset_indicator_obs] = apriori_obs
        apriori_obs_out[!subset_indicator_obs] = apriori_obs_non_assimilated

  prior_mean_out = array(0,dim=dim(x_hat))
    dimnames(prior_mean_out) = dimnames(x_hat)

  print("Done....writing inversion object output")
  return(list(posterior=list(x_hat=x_hat,P=P,inputs=inputs_list,outputs=list(modeled_obs=modeled_obs)),
  prior=list(x_hat=prior_mean_out,P=as.matrix(P_0),inputs=inputs_list,outputs=list(modeled_obs=apriori_obs_out,DFS=sum_trS,DFB=sum_trB))))

}
