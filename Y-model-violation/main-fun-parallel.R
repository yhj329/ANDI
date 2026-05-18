
library(parallel)
library(doParallel)

main_simu2_parallel<-function(nco,N=100,nt,ns,yparameter_s,yparameter_t,xi,
                     zparameter){
  Sys.setenv(OMP_NUM_THREADS=1)
  Sys.setenv(MKL_NUM_THREADS=1)
  Sys.setenv(BLAS_NUM_THREADS=1)
  hb=0.5*ns^(-0.2)
  
  set.seed(123)
  
  data<-replicate(N, data.gen2(ns=ns,nt=nt,xi=xi,yparameter_s = yparameter_s,
                               yparameter_t=yparameter_t,zparameter = zparameter), 
                  simplify = FALSE)
  
  p=ncol(data[[1]]$xs)+1
  q=ncol(data[[1]]$zs) 
  
  #--------------parallel parameter-----------
  # core number
  # Use one core by default if nco is not specified.
  if (is.null(nco)) nco <- 1
  cl <- makeCluster(nco)
  registerDoParallel(cl)
  
  # Export objects required by worker processes.
  
  clusterExport(cl, 
                varlist = c("data", "p", "q", "hb", "ns", "nt"),
                envir = environment())
  
  # Load required packages, source R functions, and compile C++ routines
  # on each worker. Thread numbers are restricted to avoid nested parallelism.
  clusterEvalQ(cl, {
    Sys.setenv(
      OMP_NUM_THREADS = "1",
      OPENBLAS_NUM_THREADS = "1",
      MKL_NUM_THREADS = "1",
      VECLIB_MAXIMUM_THREADS = "1"
    )
    Sys.setenv(OMP_NUM_THREADS=1)
    Sys.setenv(MKL_NUM_THREADS=1)
    Sys.setenv(BLAS_NUM_THREADS=1)
    library(pracma)
    library(sandwich)
    library(pROC)
    library(withr)
    library(MASS)
    library(Rcpp)
    
    source(file = "parallel_other_fun.R")
    sourceCpp("cpp/der.cpp")
    
    sourceCpp("cpp/weighted_sum.cpp")
    
    
    NULL
  })
  
  
  
  
  #-----------------main------
  results <- foreach(i = 1:N) %dopar% {
    xs=data[[i]]$xs
    zs=data[[i]]$zs
    ys=data[[i]]$ys
    
    mut=data[[i]]$mut
    psi=data[[i]]$psi
    var_nu=data[[i]]$var_nu
    
    glm_model_s<-glm(ys ~ xs+zs, family = binomial)
    
    betazshat<-glm_model_s$coefficients['zs']
    
    
    res_mle_i=glm_model_s$coefficients
    
    #------------------------------------initial--------------------------
    alpha_tilde=rootSolve::multiroot(f = function(alpha){
      estfvalue=est_full_simple(data = data[[i]],alpha = alpha,nu = c(mut,psi))
      return(estfvalue)
    },start = rep(0,3*p+q))$root
    
    
    A_total=rootSolve::gradient(f = function(alpha){
      estf_value=est_full_simple(data = data[[i]],alpha = alpha,nu = c(mut,psi))
      return(estf_value)},x =alpha_tilde )
    
    I_middle=middle_f2_simple(data = data[[i]],alpha = alpha_tilde,nu = c(mut,psi) )
    
    B_total=rootSolve::gradient(f = function(nu){
      estfvalue=est_full_simple(data = data[[i]],alpha = alpha_tilde,nu = nu)
      return(estfvalue)
    },x = c(mut,psi))
    
    sd_total=sqrt(diag(solve(A_total)%*%I_middle %*%t(solve(A_total))/nrow(xs) +
                         solve(A_total)%*%B_total%*%var_nu%*%t(B_total) %*%t(solve(A_total)))   )
    
    
    
    
    #---------------------------auc-----------------------------
    
    ####model with xz####
    p=ncol(xs)+1
    q=ncol(zs)
    #alpha_tilde=mean_alpha
    betazt=alpha_tilde[(p+1):(p+q)]
    XI=alpha_tilde[(p+q+1):(2*p+q)]
    betaxt=alpha_tilde[(2*p+q+1):(3*p+q)]
    auchat_xz_rank=auc_prob(XI = XI,beta_true = c(betaxt,betazt),beta_model =  c(betaxt,betazt),data = data[[i]])
    res_auc_xz_rank_i=auchat_xz_rank
    
    auchat_xz_kernel=auc_prob_kernel(hb = 0.5*ns^(-0.2),XI = XI,beta_true = c(betaxt,betazt),beta_model =  c(betaxt,betazt),data = data[[i]])
    res_auc_xz_kernel_i=auchat_xz_kernel
    
    ##var
    
    
    der_AUC_xz_alpha= der_auc_kernel_xz_alpha(auchat = auchat_xz_kernel,alpha = alpha_tilde,hb = 0.5*ns^(-0.2),data = data[[i]] )
    
    
    IF_source_xz=IF_kernel_source_xz(auchat = auchat_xz_kernel,der_AUC_alpha = der_AUC_xz_alpha,A_total = A_total,
                                     alpha = alpha_tilde,nu =c(mut,psi),hb = 0.5*ns^(-0.2),data = data[[i]] )
    vector_target_xz=-t(B_total) %*%t(solve(A_total))%*%drop(der_AUC_xz_alpha)
    var_target_xz= t(vector_target_xz)%*%var_nu%*%vector_target_xz
    var_source=var(IF_source_xz)
    sd_auc_xz_kernel_i=sqrt(var_source/ns +var_target_xz)
    
    
    ####model only x####
    
    betazt=alpha_tilde[(p+1):(p+q)]
    XI=alpha_tilde[(p+q+1):(2*p+q)]
    betaxt=alpha_tilde[(2*p+q+1):(3*p+q)]
    beta_true=c(betaxt,betazt)
    
    psi=psi
    beta_model=c(psi,0)
    
    auchat_x_rank=auc_prob(XI = XI,beta_true = c(betaxt,betazt),beta_model = beta_model,data = data[[i]])
    res_auc_x_rank_i=auchat_x_rank
    
    auchat_x_kernel=auc_prob_kernel(hb = 0.5*ns^(-0.2),XI = XI,beta_true = c(betaxt,betazt),beta_model = beta_model,data = data[[i]])
    res_auc_x_kernel_i=auchat_x_kernel
    
    ##var
    
    der_AUC_x_alpha=der_auc_kernel_x_alpha(auchat = auchat_x_kernel,alpha = alpha_tilde,nu = c(mut,psi),hb = 0.5*ns^(-0.2),data = data[[i]])
    der_AUC_x_nu=der_auc_kernel_x_nu(auchat = auchat_x_kernel,alpha = alpha_tilde,nu = c(mut,psi),hb = 0.5*ns^(-0.2),data=data[[i]])
    
    IF_source_x=IF_kernel_source_x(auchat = auchat_x_kernel,der_AUC_alpha = der_AUC_x_alpha,A_total = A_total,
                                   alpha = alpha_tilde,nu =c(mut,psi),hb = 0.5*ns^(-0.2),data = data[[i]] )
    
    vector_target_x=-t(B_total) %*%t(solve(A_total))%*%drop(der_AUC_x_alpha)+drop(der_AUC_x_nu)
    var_target= t(vector_target_x)%*%var_nu%*%vector_target_x
    var_source=var(IF_source_x)
    sd_auc_x_kernel_i=sqrt(var_source/ns +var_target)
    
    ####diff####
    
    res_auc_diff_rank_i=auchat_xz_rank-auchat_x_rank
    res_auc_diff_kernel_i=auchat_xz_kernel - auchat_x_kernel
    
    var_source=var(IF_source_xz- IF_source_x)
    vector_target_diff=vector_target_xz- vector_target_x
    var_target= t(vector_target_diff)%*%var_nu%*%vector_target_diff
    sd_auc_diff_kernel_i=sqrt(var_source/ns +var_target)
    
    #-----------source auc -----------------
    
    # B 200 500 1000
    # quantile and normal  
    
    
    a1=with_seed(123,boot_auc_matrix_ci(ys=ys,X=xs,Z=zs,Blist=c(200,500,1000)))
    #a1=with_seed(123,boot_auc_matrix_ci(ys=y,X=X,Z=Z,Blist=c(200,500,1000)))
    res_auc_naive_B200=a1[[1]]
    res_auc_naive_B500=a1[[2]]
    res_auc_naive_B1000=a1[[3]]
    
    
    gc()
    
    
    #-----------result-------
    out <- list(
      res_mle = res_mle_i,
      res_trans = alpha_tilde,
      sd_mle = sqrt(diag(vcov(glm_model_s))),
      sd_trans = sd_total,
      
      res_auc_xz_rank = res_auc_xz_rank_i, res_auc_xz_kernel = res_auc_xz_kernel_i,
      sd_auc_xz_kernel = sd_auc_xz_kernel_i,
      res_auc_x_rank = res_auc_x_rank_i, res_auc_x_kernel = res_auc_x_kernel_i, 
      sd_auc_x_kernel = sd_auc_x_kernel_i,
      res_auc_diff_rank = res_auc_diff_rank_i, res_auc_diff_kernel = res_auc_diff_kernel_i,
      sd_auc_diff_kernel = sd_auc_diff_kernel_i,
      
      res_auc_naive_B200 = res_auc_naive_B200, res_auc_naive_B500 = res_auc_naive_B500,
      res_auc_naive_B1000= res_auc_naive_B1000
    )
    return(out)
    
  }
  
  # stop cluster
  stopCluster(cl)
  
  #-----------summary results--------------
  res_mle <- do.call(rbind, lapply(results, function(x) x$res_mle))
  res_trans <- do.call(rbind, lapply(results, function(x) x$res_trans))
  sd_mle <- do.call(rbind, lapply(results, function(x) x$sd_mle))
  sd_trans <- do.call(rbind, lapply(results, function(x) x$sd_trans))
  
  res_auc_xz_rank <- sapply(results, function(x) x$res_auc_xz_rank)
  res_auc_xz_kernel <- sapply(results, function(x) x$res_auc_xz_kernel)
  sd_auc_xz_kernel <- sapply(results, function(x) x$sd_auc_xz_kernel)
  
  res_auc_x_rank <- sapply(results, function(x) x$res_auc_x_rank)
  res_auc_x_kernel <- sapply(results, function(x) x$res_auc_x_kernel)
  sd_auc_x_kernel <- sapply(results, function(x) x$sd_auc_x_kernel)
  
  res_auc_diff_rank <- sapply(results, function(x) x$res_auc_diff_rank)
  res_auc_diff_kernel <- sapply(results, function(x) x$res_auc_diff_kernel)
  sd_auc_diff_kernel <- sapply(results, function(x) x$sd_auc_diff_kernel)
  
  res_auc_naive_B200_est <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B200$est))
  res_auc_naive_B200_sd <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B200$sd))
  res_auc_naive_B200_ql <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B200$ql))
  res_auc_naive_B200_qu <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B200$qu))
  
  res_auc_naive_B500_est <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B500$est))
  res_auc_naive_B500_sd <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B500$sd))
  res_auc_naive_B500_ql <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B500$ql))
  res_auc_naive_B500_qu <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B500$qu))
  
  res_auc_naive_B1000_est <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B1000$est))
  res_auc_naive_B1000_sd <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B1000$sd))
  res_auc_naive_B1000_ql <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B1000$ql))
  res_auc_naive_B1000_qu <- do.call(rbind, lapply(results, function(x) x$res_auc_naive_B1000$qu))
  
  return(list(
    res_mle = res_mle, res_trans = res_trans,
    sd_mle = sd_mle, sd_trans = sd_trans,
    res_auc_xz_rank = res_auc_xz_rank, res_auc_xz_kernel = res_auc_xz_kernel,
    sd_auc_xz_kernel = sd_auc_xz_kernel,
    res_auc_x_rank = res_auc_x_rank, res_auc_x_kernel = res_auc_x_kernel,
    sd_auc_x_kernel = sd_auc_x_kernel,
    res_auc_diff_rank = res_auc_diff_rank, res_auc_diff_kernel = res_auc_diff_kernel,
    sd_auc_diff_kernel = sd_auc_diff_kernel,
    
    
    res_auc_naive_B200_est = res_auc_naive_B200_est,
    res_auc_naive_B200_sd = res_auc_naive_B200_sd,
    res_auc_naive_B200_ql = res_auc_naive_B200_ql,
    res_auc_naive_B200_qu = res_auc_naive_B200_qu,
    
    res_auc_naive_B500_est = res_auc_naive_B500_est,
    res_auc_naive_B500_sd = res_auc_naive_B500_sd,
    res_auc_naive_B500_ql = res_auc_naive_B500_ql,
    res_auc_naive_B500_qu = res_auc_naive_B500_qu,
    
    res_auc_naive_B1000_est = res_auc_naive_B1000_est,
    res_auc_naive_B1000_sd = res_auc_naive_B1000_sd,
    res_auc_naive_B1000_ql = res_auc_naive_B1000_ql,
    res_auc_naive_B1000_qu = res_auc_naive_B1000_qu,
    
    param_list=list(nt=nt,ns=ns,betat=betat,betas=betas,xi=xi,
                    zparameter=zparameter)
  ))
}



