library(pracma)
library(sandwich)
library(nloptr)
library(pROC)
library(withr)
library(MASS)

evalua<-function(res,p=4,q=1,beta_t){
  
  # targetbeta=colMeans(result$beta_cmle)
  
  # tauc[1]=mean(result$auchat_cmle[,1])
  # tauc[2]=mean(result$auchat_reduce[,1])
  # tauc[3]=mean(result$auchat_dif[,1])
  ind_t=c((2*p+q+1):(3*p+q),p+q)
  total_comp<-data.frame(c("S1","","","",""),
                         c("$\\beta_{T,X_0}$","$\\beta_{T,X_1}$","$\\beta_{T,X_2}$","$\\beta_{T,X_3}$","$\\beta_{T,Z}$"),              
                         beta_t,
                         colMeans(res$res_mle)-beta_t,
                         sqrt(diag(var(res$res_mle))),
                         colMeans(res$sd_mle),
                         colMeans(t( (t(res$res_mle)-beta_t)^2)),
                         rowMeans(abs(t(res$res_mle)-beta_t)/t(res$sd_mle)<qnorm(0.975)),
                         colMeans(res$res_trans[,ind_t])-beta_t,
                         sqrt(diag(var(res$res_trans[,ind_t]))),
                         colMeans(res$sd_trans[,ind_t]),
                         colMeans(t( (t(res$res_trans[,ind_t])-beta_t)^2)),
                         rowMeans(abs(t(res$res_trans[,ind_t])-beta_t)/t(res$sd_trans[,ind_t])<qnorm(0.975))
                         #c(rowMeans(abs(t(result$beta_cmle[,1:5])-targetbeta[1:5])/t(result$asyd_cml)<qnorm(0.975)),0,mean(abs(result$auchat_cmle[,1]-tauc[1])/sqrt(result$auchat_cmle[,2])<qnorm(0.975) ),mean(abs(result$auchat_dif[,1]-mean(result$auchat_dif[,1]))/sqrt(result$auchat_dif[,2])<qnorm(0.975) ) )
                         
  )
  auctrue=auc_true_est
  auc_xz_true_est=auc_true_est[1]
  auc_x_true_est=auc_true_est[2]
  auc_diff_true_est=auc_true_est[3]
  other_metric<-data.frame(c("","",""),
                           c("$AUC_{xz}$","$AUC_{x}$","$AUC_{diff}$"),
                           auctrue,
                           {## xz
                             rbind(
                               c(mean(res$res_auc_xz_kernel)- auc_xz_true_est,
                                 sd(res$res_auc_xz_kernel),
                                 mean(res$sd_auc_xz_kernel),
                                 mean((res$res_auc_xz_kernel- auc_xz_true_est)^2),
                                 mean(abs(res$res_auc_xz_kernel- auc_xz_true_est)/res$sd_auc_xz_kernel<1.96)),
                               ## only x
                               c(mean(res$res_auc_x_kernel)- auc_x_true_est,
                                 sd(res$res_auc_x_kernel),
                                 mean(res$sd_auc_x_kernel),
                                 mean((res$res_auc_x_kernel- auc_x_true_est)^2),
                                 mean(abs(res$res_auc_x_kernel- auc_x_true_est)/res$sd_auc_x_kernel<1.96)),
                               ## difference
                               c(mean(res$res_auc_diff_kernel)- auc_diff_true_est,
                                 sd(res$res_auc_diff_kernel),
                                 mean(res$sd_auc_diff_kernel),
                                 mean((res$res_auc_diff_kernel- auc_diff_true_est)^2),
                                 mean(abs(res$res_auc_diff_kernel- auc_diff_true_est)/res$sd_auc_diff_kernel<1.96))
                             )}
                           
  )
  #auctrue=auctrue[c(2,1,3)]
  other_metric<-cbind(other_metric,
                      #B200
                      colMeans(res$res_auc_naive_B200_est)-auctrue,
                      sqrt(diag(var(res$res_auc_naive_B200_est))),
                      colMeans(res$res_auc_naive_B200_sd),
                      colMeans(t( (t(res$res_auc_naive_B200_est)-auctrue)^2)),
                      rowMeans(abs(t(res$res_auc_naive_B200_est)-auctrue)/t(res$res_auc_naive_B200_sd)<qnorm(0.975)),
                      rowMeans(auctrue > t(res$res_auc_naive_B200_ql) &
                                 (auctrue < t(res$res_auc_naive_B200_qu))),
                      
                      #B500
                      colMeans(res$res_auc_naive_B500_est)-auctrue,
                      sqrt(diag(var(res$res_auc_naive_B500_est))),
                      colMeans(res$res_auc_naive_B500_sd),
                      colMeans(t( (t(res$res_auc_naive_B500_est)-auctrue)^2)),
                      rowMeans(abs(t(res$res_auc_naive_B500_est)-auctrue)/t(res$res_auc_naive_B500_sd)<qnorm(0.975)),
                      rowMeans(auctrue > t(res$res_auc_naive_B500_ql) &
                                 (auctrue < t(res$res_auc_naive_B500_qu))),
                      #B1000
                      colMeans(res$res_auc_naive_B1000_est)-auctrue,
                      sqrt(diag(var(res$res_auc_naive_B1000_est))),
                      colMeans(res$res_auc_naive_B1000_sd),
                      colMeans(t( (t(res$res_auc_naive_B1000_est)-auctrue)^2)),
                      rowMeans(abs(t(res$res_auc_naive_B1000_est)-auctrue)/t(res$res_auc_naive_B1000_sd)<qnorm(0.975)),
                      rowMeans(auctrue > t(res$res_auc_naive_B1000_ql) &
                                 (auctrue < t(res$res_auc_naive_B1000_qu)))
                      )
  colnames(other_metric)<-c("setting","para","True","Bias","SE","ASE","MSE","CR",
                            "Bias","SE","ASE","MSE","CR-norm","CR-quantile",
                            "Bias","SE","ASE","MSE","CR-norm","CR-quantile",
                            "Bias","SE","ASE","MSE","CR-norm","CR-quantile")
  
  
  
  
  
  colnames(total_comp)<-c("setting","para","True","Bias","SE","ASE","MSE","CR","Bias","SE","ASE","MSE","CR")
  #rownames(total_comp)<-c("$\\beta_{T,X_0}$","$\\beta_{T,X_1}$","$\\beta_{T,X_2}$","$\\beta_{T,X_3}$","$\\beta_{T,X_4}$","$\\beta_{T,Z^c}$")
  
  return(list(other_metric,total_comp))
  
}







data.gen2 <- function(ns, nt, 
                      betas, betat, 
                      mu_xt=rep(0,3), sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), xi, 
                      zparameter = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=0.7, include_interaction=FALSE),
                      threshold_p=0) {
  
  p1 <- length(mu_xt)
  
  # -------------------
  # Target: generate (X_{-p}, latent X_p) together
  # -------------------
  xt_cont <- mvrnorm(nt, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  
  xp_t <- as.numeric(xt_cont[,p1] > threshold_p)
  xt <- cbind(xt_cont[,-p1], xp_t)
  
  # Z | X
  beta_z <- zparameter$beta
  sigma_z <- zparameter$sigma
  if(zparameter$include_interaction){
    #
  } else {
    zt <- cbind(1, xt) %*% beta_z[1:(p1+1)] + rnorm(nt,0,sigma_z)
  }
  
  # Y | X, Z
  prob_yt <- 1/(1+exp(-cbind(1, xt, zt) %*% betat))
  yt <- rbinom(nt,1,prob_yt)
  
  # -------------------
  # Source: tilting for X
  # -------------------
  # adjust mean for continuous X_{-p}
  adjust_p=drop(sigma_xt[p1,-p1]) %*% xi[-p1]
  mu_source <- mu_xt - c(drop(sigma_xt[-p1,-p1] %*% xi[-p1]),adjust_p)
  
  M = max(1, exp(-xi[p1]))
  p_accept=( mean(xp_t)*exp(-xi[p1])+1-mean(xp_t) )/M
  n_temp=ceiling(1.5 * ns / p_accept)
  # 连续部分的提议分布
  x_cont_candidate  <- mvrnorm(n_temp, mu = mu_source, Sigma = sigma_xt)  # includes latent X_p
  
  # 离散部分 threshold latent X_p
  xp_candidate <- as.numeric(x_cont_candidate[,p1] >  threshold_p)
  
  
  # rejection sampling for last variable if xi[p1] != 0
  
    u <- runif(n_temp)
    accept <- u < exp(-xi[p1]*xp_candidate)/M
    x_cont_s <- x_cont_candidate[accept,,drop=FALSE]
    xp_s <- xp_candidate[accept]
  
  
  xs <- cbind(x_cont_s[,-p1,drop=FALSE], xp_s)
  xs=xs[1:ns,]
  # Z | X in target
  if(zparameter$include_interaction){
    #
  } else {
    zs <- cbind(1, xs) %*% beta_z[1:(p1+1)] + rnorm(ns,0,sigma_z)
  }
  
  # Y | X, Z in target
  prob_ys <- 1/(1+exp(-cbind(1, xs, zs) %*% betas))
  ys <- rbinom(nrow(xs),1,prob_ys)
  
  
  #-----------------summary-----------
  yt_xt_model <- glm(yt ~ xt, family = binomial)
  
  
  pred_p=drop(cbind(1,xt)%*%coef(yt_xt_model))
  pred_p=1/(1+exp(-pred_p))
  
  A_ma=t(cbind(1,xt))%*%(cbind(1,xt)*pred_p*(1-pred_p))/nrow(xt)
  
  IF_psi=(cbind(1,xt)*(yt-pred_p)) %*% t(solve(A_ma))
  var_nu=cov(cbind(1,xt,IF_psi))/nrow(xt)
  
  
  #n_ma=ncol(xt)+length(coef(yt_xt_model))
  #matrix(0,nrow=n_ma,ncol=n_ma)
  #-----------results-------------
  return(list(xs=xs,zs=zs,ys=ys,mut=colMeans(cbind(1,xt)),
              psi=coef(yt_xt_model) ,var_nu=var_nu,signal_zs=sigma_z^2/var(zs),signal_zt=sigma_z^2/var(zt)))
  
  
  #return(list(xs=xs, ys=ys, zs=zs,xt=xt, yt=yt, zt=zt))
}











data.gen2_full_ratio<- function(ns, nt, 
                      betas, betat, 
                      mu_xt=rep(0,3), sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), 
                      xi, 
                      zparameter = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=0.7, include_interaction=FALSE),
                      threshold_p=0) {
  
  p1 <- length(mu_xt)
  
  # -------------------
  # Target: generate (X_{-p}, latent X_p) together
  # -------------------
  xt_cont <- mvrnorm(nt, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  
  xp_t <- as.numeric(xt_cont[,p1] > threshold_p)
  xt <- cbind(xt_cont[,-p1], xp_t)
  
  # Z | X
  beta_z <- zparameter$beta
  sigma_z <- zparameter$sigma
  if(zparameter$include_interaction){
    #
  } else {
    zt <- cbind(1, xt) %*% beta_z[1:(p1+1)] + rnorm(nt,0,sigma_z)
  }
  
  # Y | X, Z
  prob_yt <- 1/(1+exp(-cbind(1, xt, zt) %*% betat))
  yt <- rbinom(nt,1,prob_yt)
  
  # -------------------
  # Source: tilting for X
  # -------------------
  # adjust mean for continuous X_{12}
  adjust_p=drop(sigma_xt[p1,-p1]) %*% xi[-p1]
  mu_source <- mu_xt - c(drop(sigma_xt[-p1,-p1] %*% xi[-p1]),adjust_p)
  
  M = max(1, exp(-xi[p1]))
  p_accept=( mean(xp_t)*exp(-xi[p1])+1-mean(xp_t) )/M
  n_temp=ceiling(1.5 * ns / p_accept)
  # 连续部分的提议分布
  x_cont_candidate  <- mvrnorm(n_temp, mu = mu_source, Sigma = sigma_xt)  # includes latent X_p
  
  # 离散部分 threshold latent X_p
  xp_candidate <- as.numeric(x_cont_candidate[,p1] >  threshold_p)
  
  
  # rejection sampling for last variable if xi[p1] != 0
  
  u <- runif(n_temp)
  accept <- u < exp(-xi[p1]*xp_candidate)/M
  x_cont_s <- x_cont_candidate[accept,,drop=FALSE]
  xp_s <- xp_candidate[accept]
  
  
  xs <- cbind(x_cont_s[,-p1,drop=FALSE], xp_s)
  xs=xs[1:ns,]
  # Z | X in target
  if(zparameter$include_interaction){
    #
  } else {
    zs <- cbind(1, xs) %*% beta_z[1:(p1+1)] + rnorm(ns,0,sigma_z)
  }
  
  # Y | X, Z in target
  prob_ys <- 1/(1+exp(-cbind(1, xs, zs) %*% betas))
  ys <- rbinom(nrow(xs),1,prob_ys)
  
  
  #-----------------summary-----------
  yt_xt_model <- glm(yt ~ xt, family = binomial)
  
  
  pred_p=drop(cbind(1,xt)%*%coef(yt_xt_model))
  pred_p=1/(1+exp(-pred_p))
  
  A_ma=t(cbind(1,xt))%*%(cbind(1,xt)*pred_p*(1-pred_p))/nrow(xt)
  
  IF_psi=(cbind(1,xt)*(yt-pred_p)) %*% t(solve(A_ma))
  var_nu=cov(cbind(1,xt,IF_psi))/nrow(xt)
  
  
  #n_ma=ncol(xt)+length(coef(yt_xt_model))
  #matrix(0,nrow=n_ma,ncol=n_ma)
  #-----------results-------------
  return(list(xs=xs,zs=zs,ys=ys,mut=colMeans(cbind(1,xt)),
              psi=coef(yt_xt_model) ,var_nu=var_nu,signal_zs=sigma_z^2/var(zs),signal_zt=sigma_z^2/var(zt)))
  
  
  #return(list(xs=xs, ys=ys, zs=zs,xt=xt, yt=yt, zt=zt))
}










est_full_simple<-function(data,alpha,nu){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  mut=nu[1:p]
  psi=nu[(p+1):(2*p)]
  betaxs=delta+betaxt
  
  py_xz_s=1/(1+exp(-cbind(1,xs,zs)%*%c(betaxs,betazt)) )
  A1=cbind(1,xs,zs)*as.vector((ys-py_xz_s)) 
  
  weight=drop(exp(cbind(1,xs)%*%XI) )
  
  A2=t( t(cbind(1,xs)*weight)- mut )
  
  
  py_x_psi=drop( 1/(1+exp(-cbind(1,xs)%*%psi)) )
  py_xz=drop( 1/(1+exp(-cbind(1,xs,zs)%*%c(betaxt,betazt))) )
  
  A3=cbind(1,xs)*(weight*(py_xz- py_x_psi))
  
  f_total=cbind(A1,A2,A3)
  
  return(   colMeans(f_total) )
}
middle_f2_simple<-function(data,alpha,nu){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  mut=nu[1:p]
  psi=nu[(p+1):(2*p)]
  betaxs=delta+betaxt
  
  py_xz_s=1/(1+exp(-cbind(1,xs,zs)%*%c(betaxs,betazt)) )
  A1=cbind(1,xs,zs)*as.vector((ys-py_xz_s)) 
  
  weight=drop(exp(cbind(1,xs)%*%XI) )
  
  A2=t( t(cbind(1,xs)*weight)- mut )
  
  
  py_x_psi=drop( 1/(1+exp(-cbind(1,xs)%*%psi)) )
  py_xz=drop( 1/(1+exp(-cbind(1,xs,zs)%*%c(betaxt,betazt))) )
  
  A3=cbind(1,xs)*(weight*(py_xz- py_x_psi))
  
  f_total=cbind(A1,A2,A3)
  
  return( t(f_total)%*%f_total/nrow(f_total) )
}
est_f_no_sum<-function(data,alpha,nu){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  mut=nu[1:p]
  psi=nu[(p+1):(2*p)]
  betaxs=delta+betaxt
  
  py_xz_s=1/(1+exp(-cbind(1,xs,zs)%*%c(betaxs,betazt)) )
  A1=cbind(1,xs,zs)*as.vector((ys-py_xz_s)) 
  
  weight=drop(exp(cbind(1,xs)%*%XI) )
  
  A2=t( t(cbind(1,xs)*weight)- mut )
  
  
  py_x_psi=drop( 1/(1+exp(-cbind(1,xs)%*%psi)) )
  py_xz=drop( 1/(1+exp(-cbind(1,xs,zs)%*%c(betaxt,betazt))) )
  
  A3=cbind(1,xs)*(weight*(py_xz- py_x_psi))
  
  f_total=cbind(A1,A2,A3)
  
  return( f_total )
}

#------------auc x z---------------------------------



fast_weighted_sum <- function(s, w1, w0) {
  n <- length(s)
  
  # 1. 获取排序索引
  order_idx <- order(s)  # 递增排序，order_idx 使 s[order_idx] 递增
  
  # 2. 重新排列 w1, w0
  s_sorted <- s[order_idx]
  w1_sorted <- w1[order_idx]
  w0_sorted <- w0[order_idx]
  
  # 3. 计算 w0 的前缀和
  w0_prefix_sum <- c(0, cumsum(w0_sorted))  # 长度 n+1，方便索引
  
  # 3. 计算结果
  result <- 0
  i <- 1
  while (i <= n) {
    # 找到相同 s[i] 的区间 [i, j]
    j <- i
    while (j <= n && s_sorted[j] == s_sorted[i]) j <- j + 1
    j <- j - 1  # j 是最后一个等于 s[i] 的索引
    
    # 处理所有 s[i] > s[j] 的部分
    if (i > 1) {
      result <- result + sum(w1_sorted[i:j]) * w0_prefix_sum[i]
    }
    
    # 处理所有 s[i] = s[j] 的部分 (二项式求和)
    equal_w1_sum <- sum(w1_sorted[i:j])
    equal_w0_sum <- sum(w0_sorted[i:j])
    result <- result + 0.5 * equal_w1_sum * equal_w0_sum
    
    # 移动到下一个不同的 s[i]
    i <- j + 1
  }
  
  return(result)
}

auc_prob<-function(XI,beta_true,beta_model,data){
  
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  weight=drop(exp(cbind(1,xs)%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  
  sum_auc<-fast_weighted_sum(s =score,w1 = pi1,w0 = pi0 )
  #compute_sum(out_ma_model=out_ma_model, w0=w0, w1=w1) 
  
  sum_auc=(sum_auc )/(sum(pi1)*sum(pi0)) #/P(Y=1)*P(Y=0) 
  return(sum_auc)
}


weighted_sum_kernel <- function(hb,s, w1, w0) {
  n <- length(s)
  
  # 构造得分差矩阵 (n1 x n0)
  diff_mat <- outer(s, s, "-") / hb
  
  # 计算Phi值矩阵
  Phi_mat <- pnorm(diff_mat)
  
  # 权重矩阵 = w_pos[i] * w_neg[j]
  w_mat <- outer(w1, w0, "*")
  
  # AUC = sum_{i,j} w1[i] w0[j] Phi(...)
  auc_est <- sum(w_mat * Phi_mat)
  
  return(auc_est)
}

auc_prob_kernel<-function(hb,XI,beta_true,beta_model,data){
  
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  weight=drop(exp(cbind(1,xs)%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  
  sum_auc<-weighted_sum_kernel(hb = hb,s =score,w1 = pi1,w0 = pi0 )
  #compute_sum(out_ma_model=out_ma_model, w0=w0, w1=w1) 
  
  sum_auc=(sum_auc )/(sum(pi1)*sum(pi0)) #/P(Y=1)*P(Y=0) 
  return(sum_auc)
}

library(Rcpp)
sourceCpp("cpp/der.cpp")

sourceCpp("cpp/weighted_sum.cpp")
## var kernel source部分
IF_kernel_source_xz<-function(auchat,der_AUC_alpha,A_total,alpha,nu,
                              hb,data ){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  
  beta_model=c(betaxt,betazt)
  beta_true=c(betaxt,betazt)
  
  weight=drop(exp(cbind(1,xs)%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  
  pi1_mean= mean(pi1)
  pi0_mean=mean(pi0)
  
  IF_1=weighted_IF_cpp(S =  score,pi1 = pi1,pi0 = pi0,hb = hb)/(pi1_mean*pi0_mean)-
    auchat*(pi0/pi0_mean+pi1/pi1_mean)
  {
    # # 构造得分差矩阵 (n1 x n0)
    # diff_mat <- outer(score, score, "-") / hb
    # 
    # # 计算Phi值矩阵
    # Phi_mat <- pnorm(diff_mat)
    # 
    # # 权重矩阵 = w_pos[i] * w_neg[j]
    # w_mat <- outer(pi1, pi0, "*")
    # IF_1=drop(colMeans(w_mat * Phi_mat+t(w_mat * Phi_mat)))/(mean(pi0)*mean(pi1))- #part1
    #   auchat*(pi0/mean(pi0)+pi1/mean(pi1)) #part2
    
    }
  IF_1=scale( IF_1,center=T,scale=F)
  
  #return(var(IF_1))
  IF_2=IF_1-drop(est_f_no_sum(data = data,alpha =alpha ,nu=nu)%*%t(solve(A_total))%*%drop(der_AUC_alpha))#part 3_1
  
  return((IF_2))
  #return(mean(IF_2^2))
}

der_auc_kernel_xz_alpha<-function(auchat,alpha,hb,data){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  
  beta_model=c(betaxt,betazt)
  beta_true=c(betaxt,betazt)
  
  
  ma_X=cbind(1,xs)
  
  weight=drop(exp(ma_X%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  pi1_mean= mean(pi1)
  pi0_mean=mean(pi0)
  
  der_XI=der_h_weighting_XI(pi1 = pi1,pi0 = pi0,S = score,x = ma_X,hb = hb)/
    (pi1_mean*pi0_mean)
  der_XI=der_XI-auchat* drop(colMeans(  ma_X*pi0/pi0_mean+ma_X*pi1/pi1_mean))
  
  der_beta=der_h_model_beta(pi1 = pi1,pi0 = pi0,S = score,py = y_p,Predictor = xzdata,hb = hb)+
    der_h_weighting_beta(pi1 = pi1,pi0 = pi0,S = score,py = y_p,W = xzdata,hb = hb)
  der_beta=der_beta/(pi1_mean*pi0_mean)
  der_beta=der_beta-auchat* drop(colMeans(  xzdata*(1-y_p)*pi1/pi1_mean-xzdata*y_p*pi0/pi0_mean))
  
  der=rep(0,length(alpha))
  der[(p+q+1):(2*p+q)]=  der_XI
  der[c((2*p+q+1):(3*p+q),(p+1):(p+q))]=der_beta
  return(der)
}
##------------------------ auc only x model----------------------------


###方差
## var kernel source部分
IF_kernel_source_x<-function(auchat,der_AUC_alpha,A_total,alpha,nu,
                             hb,data ){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  
  beta_true=c(betaxt,betazt)
  
  psi=nu[(p+1):(2*p)]
  beta_model=c(psi,0)
  
  
  weight=drop(exp(cbind(1,xs)%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  
  pi1_mean= mean(pi1)
  pi0_mean=mean(pi0)
  
  IF_1=weighted_IF_cpp(S = score,pi1 = pi1,pi0 = pi0,hb = hb)/(pi1_mean*pi0_mean)-
    auchat*(pi0/pi0_mean+pi1/pi1_mean)
  
  {
    # # 构造得分差矩阵 (n1 x n0)
    # diff_mat <- outer(score, score, "-") / hb
    # 
    # # 计算Phi值矩阵
    # Phi_mat <- pnorm(diff_mat)
    # 
    # # 权重矩阵 = w_pos[i] * w_neg[j]
    # w_mat <- outer(pi1, pi0, "*")
    # IF_1=drop(colMeans(w_mat * Phi_mat+t(w_mat * Phi_mat)))/(mean(pi0)*mean(pi1))- #part1
    #   auchat*(pi0/mean(pi0)+pi1/mean(pi1)) #part2
    }
  IF_1=scale( IF_1,center=T,scale=F)
  #return(var(IF_1))
  IF_2=IF_1-drop(est_f_no_sum(data = data,alpha =alpha ,nu=nu)%*%t(solve(A_total))%*%drop(der_AUC_alpha))#part 3_1
  
  return((IF_2))
  #return(mean(IF_2^2))
}

der_auc_kernel_x_alpha<-function(auchat,alpha,nu,hb,data){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  
  psi=nu[(p+1):(2*p)]
  beta_model=c(psi,0)
  beta_true=c(betaxt,betazt)
  
  
  ma_X=cbind(1,xs)
  
  weight=drop(exp(ma_X%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  pi1_mean= mean(pi1)
  pi0_mean=mean(pi0)
  der_XI=der_h_weighting_XI(pi1 = pi1,pi0 = pi0,S = score,x = ma_X,hb = hb)/
    (pi1_mean*pi0_mean)
  der_XI=der_XI-auchat* drop(colMeans(  ma_X*pi0/pi0_mean+ma_X*pi1/pi1_mean))
  
  der_beta=der_h_weighting_beta(pi1 = pi1,pi0 = pi0,S = score,py = y_p,W = xzdata,hb = hb)
  der_beta=der_beta/(pi1_mean*pi0_mean)
  der_beta=der_beta-auchat* drop(colMeans(  xzdata*(1-y_p)*pi1/pi1_mean-xzdata*y_p*pi0/pi0_mean))
  
  der=rep(0,length(alpha))
  der[(p+q+1):(2*p+q)]=  der_XI
  der[c((2*p+q+1):(3*p+q),(p+1):(p+q))]=der_beta
  return(der)
}


der_auc_kernel_x_nu<-function(auchat,alpha,nu,hb,data ){
  xs=data$xs
  zs=data$zs
  ys=data$ys
  p=ncol(xs)+1
  q=ncol(zs)
  delta=alpha[1:p]
  betazt=alpha[(p+1):(p+q)]
  XI=alpha[(p+q+1):(2*p+q)]
  betaxt=alpha[(2*p+q+1):(3*p+q)]
  
  psi=nu[(p+1):(2*p)]
  beta_model=c(psi,0)
  beta_true=c(betaxt,betazt)
  
  
  ma_X=cbind(1,xs)
  
  weight=drop(exp(ma_X%*%XI) )
  xzdata=cbind(1,xs,zs)
  score=drop(xzdata%*%beta_model)
  y_p=drop(xzdata%*%beta_true)
  y_p=1/(1+exp(-y_p))
  pi1=weight*y_p
  pi0=weight*(1-y_p)
  pi1_mean= mean(pi1)
  pi0_mean=mean(pi0)
  
  der_psi=der_h_model_beta(pi1 = pi1,pi0 = pi0,S = score,py = y_p,Predictor = ma_X,hb = hb)
  der_psi=der_psi/(pi1_mean*pi0_mean)
  der=rep(0,length(nu))
  der[(p+1):(2*p)]=  der_psi
  return(der)
}



#-------auc usual---------------


auc_multi<-function(tlabel,pred_list){
  if (is.null(names(pred_list))) {
    stop("pred_list 必须是一个具名 list，例如 list(M1 = pred1, M2 = pred2)")
  }
  
  # 确保tlabel是数值向量
  if (!is.numeric(tlabel)) {
    stop("tlabel 必须是数值向量")
  }
  
  
  
  
  
  # 确保pred_list中所有元素都是数值向量且长度与tlabel一致
  for (i in seq_along(pred_list)) {
    if (!is.numeric(pred_list[[i]])) {
      stop(paste0("pred_list 中的元素 '", names(pred_list)[i], "' 必须是数值向量"))
    }
    if (length(pred_list[[i]]) != length(tlabel)) {
      stop(paste0("pred_list 中的元素 '", names(pred_list)[i], "' 的长度必须与 tlabel 一致"))
    }
  }
  
  
  result_list <- list()
  
  
  # 对每个模型和每个tau计算TPR
  for (model_name in names(pred_list)) {
    predictions <- pred_list[[model_name]]
    
    auc_value <- as.numeric(pROC::auc(tlabel, predictions,quiet = TRUE))
    
    # 将结果存入列表，使用模型名作为元素名
    result_list[[model_name]] <- auc_value
  }
  
  
  
  return(result_list)
  
  
}



tpr_multi <- function(tlabel, pred_list, tau) {
  if (is.null(names(pred_list))) {
    stop("pred_list 必须是一个具名 list，例如 list(M1 = pred1, M2 = pred2)")
  }
  
  # 确保tlabel是数值向量
  if (!is.numeric(tlabel)) {
    stop("tlabel 必须是数值向量")
  }
  
  # 确保pred_list中所有元素都是数值向量且长度与tlabel一致
  for (i in seq_along(pred_list)) {
    if (!is.numeric(pred_list[[i]])) {
      stop(paste0("pred_list 中的元素 '", names(pred_list)[i], "' 必须是数值向量"))
    }
    if (length(pred_list[[i]]) != length(tlabel)) {
      stop(paste0("pred_list 中的元素 '", names(pred_list)[i], "' 的长度必须与 tlabel 一致"))
    }
  }
  
  # 确保tau是数值标量或向量
  if (!is.numeric(tau)) {
    stop("tau 必须是数值标量或向量")
  }
  
  
  result_list <- list()
  
  # 对每个模型计算TPR
  for (model_name in names(pred_list)) {
    predictions <- pred_list[[model_name]]
    
    # 应用阈值得到预测标签
    pred_labels <- as.numeric(predictions >= tau)
    
    # 计算TPR
    tp <- sum(tlabel == 1 & pred_labels == 1)  # 真阳性
    fn <- sum(tlabel == 1 & pred_labels == 0)  # 假阴性
    
    # 避免除以零
    if ((tp + fn) == 0) {
      tpr_value <- 0
    } else {
      tpr_value <- tp / (tp + fn)
    }
    
    # 将结果存入列表，使用模型名作为元素名
    result_list[[model_name]] <- tpr_value
  }
  
  return(result_list)
  
}


light_auc<- function(y, score) {
  r <- rank(score)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
#----------------------auc-boot----------------
boot_auc_matrix <- function(ys,X, B=100){
  
  n <- nrow(X)
  
  # 将数据组合为数据框
  train_data <- data.frame(ys = ys, X = X)
  
  # apparent AUC
  fit <- glm(ys ~ ., data = train_data, family = binomial(link = "logit"))
  
  pred <- predict(fit, type="response")
  auc0 <- pROC::roc(ys, pred, quiet=TRUE)$auc
  
  optimism_vec <- numeric(B)
  
  for(b in 1:B){
    idx <- sample.int(n, replace=TRUE)
    
    # bootstrap sample model
    fit_b <- glm(ys ~ .,data=train_data[idx,], family=binomial())
    
    # apparent AUC on bootstrap sample
    pred_b_app <- predict(fit_b, type="response")
    auc_b_app <- pROC::roc(train_data[idx,1], pred_b_app, quiet=TRUE)$auc
    
    # test AUC on original data
    pred_b_test <- predict(fit_b, newdata=train_data[,-1], type="response")
    
    auc_b_test <- pROC::roc(ys, pred_b_test, quiet=TRUE)$auc
    
    optimism_vec[b] <- auc_b_app - auc_b_test
  }
  
  optimism <- mean(optimism_vec)
  auc_corrected <- auc0 - optimism
  
  list(
    auc_apparent = auc0,
    optimism = optimism,
    auc_corrected = auc_corrected
  )
}


boot_auc_matrix_ci <- function(ys,X,Z, Blist=c(100)){
  
  n <- nrow(X)
  
  #-----------将数据组合为数据框--------------
  train_data <- data.frame(ys = ys, X = X,Z=Z)
  col_X<-paste0("X", 1:ncol(X))
  
  if(is.matrix(Z)) {
    col_Z <- paste0("Z", 1:ncol(Z))
  }else{col_Z <- paste0("Z", 1)}
  colnames(train_data)<-c("ys",col_X,col_Z)
  
  #--------------创建公式----------------------
  # 模型1：ys 对 X 的逻辑回归
  formula_ys_x <- as.formula(paste("ys ~", paste(col_X, collapse = " + ")))
  
  # 模型2：ys 对 X 和 Z 的逻辑回归
  formula_ys_xz <- as.formula(paste("ys ~", paste(c(col_X, col_Z), collapse = " + ")))
  
  #------------apparent AUC--------------
  # 模型1：ys ~ X
  model_ys_x <- glm(formula = formula_ys_x,
                    data = train_data,
                    family = binomial(link = "logit"))


  # 模型2：ys ~ X + Z
  model_ys_xz <- glm(formula = formula_ys_xz,
                     data = train_data,
                     family = binomial(link = "logit"))

  pred_x <- predict(model_ys_x, type="response")

  auc0 <-light_auc(y = ys,score = pred_x)
    #pROC::roc(ys, pred_x, quiet=TRUE)$auc

  pred_xz <- predict(model_ys_xz, type="response")
  auc1 <- light_auc(y = ys,score = pred_xz)
  #pROC::roc(ys, pred_xz, quiet=TRUE)$auc
  
  
  #-------------boot------------------
  Bmax=max(Blist)
  auc0_b_app<-rep(0,Bmax)
  auc0_b_test<-rep(0,Bmax)
  auc1_b_app<-rep(0,Bmax)
  auc1_b_test<-rep(0,Bmax)
  
  for(b in 1:Bmax){
    idx <- sample.int(n, replace=TRUE)

    # bootstrap sample model
    # 模型1：ys ~ X

    model_ys_x_b <- glm(formula = formula_ys_x,
                      data = train_data[idx,],
                      family = binomial(link = "logit"))


    # 模型2：ys ~ X + Z
    model_ys_xz_b <- glm(formula = formula_ys_xz,
                       data = train_data[idx,],
                       family = binomial(link = "logit"))



    # apparent AUC on bootstrap sample
    pred_x_b_app <- predict(model_ys_x_b, type="response")
    auc0_b_app[b] <- light_auc(y = ys[idx],score = pred_x_b_app)
      #pROC::roc(ys[idx], pred_x_b_app, quiet=TRUE)$auc

    pred_xz_b_app <- predict(model_ys_xz_b, type="response")
    auc1_b_app[b] <-  light_auc(y = ys[idx],score = pred_xz_b_app)
      #pROC::roc(ys[idx], pred_xz_b_app, quiet=TRUE)$auc


    # test AUC on original data
    pred_x_b_test <- predict(model_ys_x_b, newdata=train_data,type="response")
    auc0_b_test[b] <-  light_auc(y = ys,score = pred_x_b_test)
      #pROC::roc(ys, pred_x_b_test, quiet=TRUE)$auc

    pred_xz_b_test <- predict(model_ys_xz_b, newdata=train_data,type="response")
    auc1_b_test[b] <- light_auc(y = ys,score = pred_xz_b_test)
      #pROC::roc(ys, pred_xz_b_test, quiet=TRUE)$auc

  }
  
  #------------------- 储存结果-------------------
  boot_results=list()

  for(i in 1:length(Blist)){
    B=Blist[i]
    optimism0 <- mean(auc0_b_app[1:B] - auc0_b_test[1:B])
    optimism1 <- mean(auc1_b_app[1:B] - auc1_b_test[1:B])
    auc0_corrected <- auc0 - optimism0
    auc1_corrected <- auc1 - optimism1
    # naive 点估计
    est_naive=c(auc1_corrected,auc0_corrected,auc1_corrected-auc0_corrected)
    # se
    sd_naive=c(sd(auc1_b_app[1:B]),
      sd(auc0_b_app[1:B]),
      sd(auc1_b_app[1:B] - auc0_b_app[1:B]))
    # ql
    ql=c(quantile(auc1_b_app[1:B], probs = 0.025),
         quantile(auc0_b_app[1:B], probs = 0.025),
         quantile(auc1_b_app[1:B]-auc0_b_app[1:B], probs = 0.025))-c(optimism1,optimism0,optimism1-optimism0)
    # qu
    qu=c(quantile(auc1_b_app[1:B], probs = 0.975),
         quantile(auc0_b_app[1:B], probs = 0.975),
         quantile(auc1_b_app[1:B]-auc0_b_app[1:B], probs = 0.975))-c(optimism1,optimism0,optimism1-optimism0)
  
    boot_results[[i]]=list(B=B,est=est_naive,sd=sd_naive,ql=ql,qu=qu)
  }
  
  return(boot_results)
}



