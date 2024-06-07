#################################################################################################
#-- PART III:    Basic inversion wrapper
#################################################################################################

invert_clean = function(H,R_diagonal,P_0,y,H_bgd,subset_indicator_obs=NULL, 
                        DOF=FALSE, output_Kalman_Gain=FALSE, output_Infl_Matrix=FALSE,force=FALSE,
                       state_vector_true=state_vector_true)
{ 
  # H=jacob;R_diagonal=R_diagonal_in;P_0=sigma;y=y_in;H_bgd=jacob_bgd

###########################################
#-- DATA MANAGEMENT, BOOKKEEPING CODE BELOW
###########################################
    
  #-- capture full sensitivity matrix info before subsetting
  if(is.null(subset_indicator_obs)){subset_indicator_obs=rep(TRUE,length(R_diagonal))}


  #-- check for too large of optimization as it will be die on GHGHub (30GB limit)
  if(!force){
     if(sum(subset_indicator_obs) >= 0.5*(dim(H)[1])) {stop(paste("Please limit observations to",
                                                            0.5*(dim(H)[1]),"observations"))}
      }
    
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

###########################################
#-- DATA MANAGEMENT, BOOKKEEPING CODE ABOVE
###########################################
    
################################
#-- CORE INVERSION CODE BELOW
################################

  R_diagonal_inv = R_diagonal^-1
  R2_diagonal_inv = R_diagonal^-2
  
  HR = H * R_diagonal_inv
  
  print("...cross product")
  
  tH_H = crossprod(HR)
  
  rm(HR)
  gc()
  
  computation_1 = solve(P_0) + tH_H
  
  rm(tH_H)
  
  print(".. .deriving posterior covariance matrix of state, P")
  P = solve(computation_1)
    
  apriori_obs = H %*% rep(1+0,dim(H)[2])
                        #generate_observations(H=H,H_bgd=H_bgd,
                        #state_vector=rep(0,(dim(H)[2])),err_obs=NULL)

  if(sum(!assim_T_F) == 1){
      apriori_obs_non_assimilated = sum(H_non_assimilated * rep(1+0,length(H_non_assimilated))) 
    }else{
      apriori_obs_non_assimilated = H_non_assimilated %*% rep(1+0,dim(H_non_assimilated)[2]) 
    }
    
  computation_3 = matrix((apriori_obs)-y,ncol=1)
  
  HR2 = H * R2_diagonal_inv
  
  HR2 = t(HR2)

  computation_4 = HR2  %*% computation_3 #+ xtra
    
  print("...deriving posterior mean state, X_hat")
  x_hat =  - P %*% computation_4

################################
#-- CORE INVERSION CODE ABOVE
################################

################################
#-- DIAGNOSTICS BELOW
################################
  
  if(output_Kalman_Gain){
    prod1 = P_0 %*% t(H) 
    sum1 = t(prod1) * (R2_diagonal_inv)
    sum1 = t(sum1)
    sum2 = - prod1 %*% t(HR2) %*% P %*% HR2
    K = sum1 + sum2
    rm(prod1);rm(sum1);rm(sum2);
  }else{K = NULL}
  
  rm(HR2)
  
  if(output_Infl_Matrix){
    P_tH = P %*% t(H)
    infl_matrix = colSums(t(H) * P_tH) * R2_diagonal_inv
    rm(P_tH)
  }
    
  #---  Digression to calc S, the influence/sensivity matrix, DOF, and chi squares

  if(DOF){
    S3 = (R2_diagonal_inv * (H) ) * t(P%*%t(H))
    trS = apply(S3,1,sum)
    sum_trS = sum(trS)
    sum_trB = dim(P)[1] - sum(trS)
    print(paste("DOF Signal:",sum_trS))
    print(paste("DOF Background:",sum_trB))     
    }else{sum_trS=NULL;sum_trB=NULL}

  #-------------------------------------------------------------
    

  
  #-- Clean up objects
  #rm(HR2)
  rm(computation_4)
  rm(computation_3)

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

#-- End modeled obs construction
    
#-- Print chi square on posterior resids vs assumed MDM
  ch_ret = varTest(as.numeric((modeled_obs[subset_indicator_obs]-y)/R_diagonal))
  print("")
  print("************************************************************")
  print("--Chi sq test on posterior residuals relative to S_z--")
  print("************************************************************")
  print(paste("var est=",floor(ch_ret$estimate*1e3)*1e-3," CI (stand variance, chi sq test): ",
              "(",floor(ch_ret$conf.int[["LCL"]]*1e3)*1e-3,",",ceiling(ch_ret$conf.int[["UCL"]]*1e3)*1e-3,")"))

#-- Print chi square on posterior state vs truth, relative to posterior state cov
  print("")
  print("************************************************************")
  print("--Chi sq test on posterior vs truth, relative to S_xpost--")
  print("************************************************************" ) 
  cstat = t(x_hat-state_vector_true) %*% solve(P) %*% (x_hat-state_vector_true) * 1/length(x_hat)
  print(paste("chi sq stat:",cstat))
  
#-- Print chi square on posterior state vs truth, relative to posterior state cov
#  print("")
#  print("*****************************************************************")
#  print("--Chi sq test: test variance of x_hat-x_prior, relative to S_0--")
#  print("*****************************************************************")  
#  cstat = t(x_hat) %*% solve(P_0) %*% x_hat * 1/length(x_hat)
#  print(paste("chi sq stat:",cstat))

  #-- TEST chi square on posterior state vs truth, relative to posterior state cov
  print("")
  print("*****************************************************************")
  print("--Chi sq test: test variance of x_prior, relative to S_0--")
  print("*****************************************************************")  
  cstat = t(state_vector_true) %*% solve(P_0) %*% state_vector_true * 1/length(state_vector_true)
  print(paste("chi sq stat:",cstat))
  
#--  Write out the inversion "object"
  print("Done....writing inversion object output")  
  return(list(posterior=list(x_hat=x_hat,P=as.matrix(P),inputs=inputs_list,outputs=list(modeled_obs=modeled_obs)),
              prior=list(x_hat=prior_mean_out,P=as.matrix(P_0),inputs=inputs_list,outputs=list(modeled_obs=apriori_obs_out)),
              diags=list(KGAIN=K,DFS=sum_trS,DFB=sum_trB)))
             
}



invert_clean_notation = function(H,Sz_diagonal,Sx,z,H_bgd,subset_indicator_obs=NULL, 
                        DOF=FALSE, output_Kalman_Gain=FALSE, output_Infl_Matrix=FALSE,force=FALSE,
                       state_vector_true=state_vector_true)
{ 
  # H=jacob;R_diagonal=R_diagonal_in;P_0=sigma;y=y_in;H_bgd=jacob_bgd

###########################################
#-- DATA MANAGEMENT, BOOKKEEPING CODE BELOW
###########################################
    
  #-- capture full sensitivity matrix info before subsetting
  if(is.null(subset_indicator_obs)){subset_indicator_obs=rep(TRUE,length(Sz_diagonal))}


  #-- check for too large of optimization as it will be die on GHGHub (30GB limit)
  if(!force){
     if(sum(subset_indicator_obs) >= 0.5*(dim(H)[1])) {stop(paste("Please limit observations to",
                                                            0.5*(dim(H)[1]),"observations"))}
      }
    
  obs_id = dimnames(H)[[1]]
  obs_len_init = dim(H)[1]
  assim_T_F=rep(0,dim(H)[1])
  assim_T_F[subset_indicator_obs] = 1
  state_names = dimnames(H)[[2]]
  
  inputs_list = list(obs_errors=Sz_diagonal,assim_T_F=assim_T_F,obs_id=obs_id,obs=z)
  

  #-- Initial subsetting of sensitivity matrix and observations
  #-- Note sensitivity matrix H is now *subset*
    H_non_assimilated=H[!subset_indicator_obs,]
    H_bgd_non_assimilated=H_bgd[!subset_indicator_obs,]    
    H=H[subset_indicator_obs,]
    H_bgd=H_bgd[subset_indicator_obs,]
    Sz_diagonal = Sz_diagonal[subset_indicator_obs]
    z = z[subset_indicator_obs]
  #}

###########################################
#-- DATA MANAGEMENT, BOOKKEEPING CODE ABOVE
###########################################
    
################################
#-- CORE INVERSION CODE BELOW
################################

  Sz_diagonal_inv = Sz_diagonal^-1
  Sz2_diagonal_inv = Sz_diagonal^-2
  
  HSz = H * Sz_diagonal_inv
  
  print("...cross product")
  
  tH_H = crossprod(HSz)
  
  rm(HSz)
  gc()
  
  computation_1 = solve(Sx) + tH_H
  
  rm(tH_H)
  
  print(".. .deriving posterior covariance matrix of state, Sx_post")
  Sx_post = solve(computation_1)
    
  apriori_obs = H %*% rep(1+0,dim(H)[2])

  if(sum(!assim_T_F) == 1){
      apriori_obs_non_assimilated = sum(H_non_assimilated * rep(1+0,length(H_non_assimilated))) 
    }else{
      apriori_obs_non_assimilated = H_non_assimilated %*% rep(1+0,dim(H_non_assimilated)[2]) 
    }
    
  computation_3 = matrix((apriori_obs)-z,ncol=1)
  
  HSz2 = H * Sz2_diagonal_inv
  
  HSz2 = t(HSz2)

  computation_4 = HSz2  %*% computation_3 #+ xtra
    
  print("...deriving posterior mean state, X_hat")
  x_hat =  - Sx_post %*% computation_4

################################
#-- CORE INVERSION CODE ABOVE
################################

################################
#-- DIAGNOSTICS BELOW
################################
  
  if(output_Kalman_Gain){
    prod1 = Sx %*% t(H) 
    sum1 = t(prod1) * (Sz2_diagonal_inv)
    sum1 = t(sum1)
    sum2 = - prod1 %*% t(HSz2) %*% Sx_post %*% HSz2
    K = sum1 + sum2
    rm(prod1);rm(sum1);rm(sum2);
  }else{K = NULL}
  
  rm(HSz2)
  
  if(output_Infl_Matrix){
    Sx_post_tH = Sx_post %*% t(H)
    infl_matrix = colSums(t(H) * Sx_post_tH) * Sz2_diagonal_inv
    rm(Sx_post_tH)
  }
    
  #---  Digression to calc S, the influence/sensivity matrix, DOF, and chi squares

  if(DOF){
    S3 = (Sz2_diagonal_inv * (H) ) * t(Sx_post%*%t(H))
    trS = apply(S3,1,sum)
    sum_trS = sum(trS)
    sum_trB = dim(Sx_post)[1] - sum(trS)
    print(paste("DOF Signal:",sum_trS))
    print(paste("DOF Background:",sum_trB))     
    }else{sum_trS=NULL;sum_trB=NULL}

  #-------------------------------------------------------------
    

  
  #-- Clean up objects
  rm(computation_4)
  rm(computation_3)

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

#-- End modeled obs construction
    
#-- Print chi square on posterior resids vs assumed MDM
  ch_ret = varTest(as.numeric((modeled_obs[subset_indicator_obs]-z)/Sz_diagonal))
  print("")
  print("************************************************************")
  print("--Chi sq test on posterior residuals relative to S_z--")
  print("************************************************************")
  print(paste("var est=",floor(ch_ret$estimate*1e3)*1e-3," CI (stand variance, chi sq test): ",
              "(",floor(ch_ret$conf.int[["LCL"]]*1e3)*1e-3,",",ceiling(ch_ret$conf.int[["UCL"]]*1e3)*1e-3,")"))

#-- Print chi square on posterior state vs truth, relative to posterior state cov
  print("")
  print("************************************************************")
  print("--Chi sq test on posterior vs truth, relative to S_xpost--")
  print("************************************************************" ) 
  cstat = t(x_hat-state_vector_true) %*% solve(Sx_post) %*% (x_hat-state_vector_true) * 1/length(x_hat)
  print(paste("chi sq stat:",cstat))

  #-- TEST chi square on posterior state vs truth, relative to posterior state cov
  print("")
  print("*****************************************************************")
  print("--Chi sq test: test variance of x_prior, relative to S_0--")
  print("*****************************************************************")  
  cstat = t(state_vector_true) %*% solve(Sx) %*% state_vector_true * 1/length(state_vector_true)
  print(paste("chi sq stat:",cstat))
  
#--  Write out the inversion "object"
  print("Done....writing inversion object output")  
  return(list(posterior=list(x_hat=x_hat,Sx_post=as.matrix(Sx_post),inputs=inputs_list,outputs=list(modeled_obs=modeled_obs)),
              prior=list(x_hat=prior_mean_out,Sx=as.matrix(Sx),inputs=inputs_list,outputs=list(modeled_obs=apriori_obs_out)),
              diags=list(KGAIN=K,DFS=sum_trS,DFB=sum_trB)))
             
}