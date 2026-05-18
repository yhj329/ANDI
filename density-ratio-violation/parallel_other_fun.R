library(pracma)
library(sandwich)
library(nloptr)
library(pROC)
library(withr)
library(MASS)



#---------------density ratio misspecification data generation--------------
data.gen2 <- function(ns, nt,betas, betat,
                      mu_xt=rep(0,3), sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), xi,
                      zparameter,threshold_p=0) {
  
  
  
  # -------------------
  # Target: generate (X_{-p}, latent X_p) together
  # -------------------
  xt_cont <- mvrnorm(nt, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  
  xp_t <- as.numeric(xt_cont[,3] > threshold_p)
  xt <- cbind(xt_cont[,-3], xp_t)
  
  # Z | X
  beta_z <- zparameter$beta
  sigma_z <- zparameter$sigma_z
  zt <- cbind(1, xt) %*% beta_z[1:(3+1)] + rnorm(nt,0,sigma_z)
  
  
  # Y | X, Z
  prob_yt <- 1/(1+exp(-cbind(1, xt, zt) %*% betat))
  yt <- rbinom(nt,1,prob_yt)
  
  # -------------------
  # Source: tilting for X
  # -------------------
  
  #------------tilting parameter--------------------------
  xi_main=xi$xi_main
  xi_x_cross=xi$xi_x_cross
  xi_zx_cross=xi$xi_zx_cross
  
  xi_prop=c(0,xi_main,xi_x_cross,xi_zx_cross)
  matr_ratio_full_term_t=cbind(1,xt,zt,xt[,1]*xt[,2],xt[,1]*xt[,3],xt[,2]*xt[,3],zt*xt[,1],zt*xt[,2],zt*xt[,3])
  
  seq_ratio_t=exp(-(matr_ratio_full_term_t%*%xi_prop) )
  M=1.5*quantile(seq_ratio_t,0.995)
  
  n_prop=floor(1.5*ns/mean(seq_ratio_t/M)) 

  #------------------------proposal dist (target)------------------------
  xprop_cont <- mvrnorm(n_prop, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  
  x3_prop <- as.numeric(xprop_cont[,3] > threshold_p)
  x_prop <- cbind(xprop_cont[,-3], x3_prop)
  
  # Z | X in target
  z_prop <- cbind(1, x_prop) %*% (beta_z[1:4]) + rnorm(n_prop,0,sigma_z)
  
  
  
  matr_ratio_full_term_prop=cbind(1,x_prop,z_prop,x_prop[,1]*x_prop[,2],x_prop[,1]*x_prop[,3],x_prop[,2]*x_prop[,3],z_prop*x_prop[,1],z_prop*x_prop[,2],z_prop*x_prop[,3])
  

  
  seq_ratio=exp(-(matr_ratio_full_term_prop%*%xi_prop) )
  
  #---------------rejection sampling--------------------
  
  u <- runif(n_prop)
  accept <- u < (seq_ratio/M)
  xs <- x_prop[accept,,drop=FALSE]
  zs <- z_prop[accept]
  
  
  xs=xs[1:ns,]
  zs=zs[1:ns]
  
  gc()
  
  # Y | X, Z 
  prob_ys <- 1/(1+exp(-cbind(1, xs, zs) %*% betas))
  ys <- rbinom(nrow(xs),1,prob_ys)
  
  
  #-----------------summary-----------
  yt_xt_model <- glm(yt ~ xt, family = binomial)
  
  
  pred_p=drop(cbind(1,xt)%*%coef(yt_xt_model))
  pred_p=1/(1+exp(-pred_p))
  
  A_ma=t(cbind(1,xt))%*%(cbind(1,xt)*pred_p*(1-pred_p))/nrow(xt)
  
  IF_psi=(cbind(1,xt)*(yt-pred_p)) %*% t(solve(A_ma))
  var_nu=cov(cbind(1,xt,IF_psi))/nrow(xt)
  
  

  #-----------results-------------
  return(list(xs=xs,zs=as.matrix(zs,ncol=1),ys=ys,mut=colMeans(cbind(1,xt)),
              psi=coef(yt_xt_model) ,var_nu=var_nu,signal_zs=sigma_z^2/var(zs),signal_zt=sigma_z^2/var(zt)))
  
  
}
#-------------------estimating equation-------------------------------
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
  
  
  order_idx <- order(s)  
  
  
  s_sorted <- s[order_idx]
  w1_sorted <- w1[order_idx]
  w0_sorted <- w0[order_idx]
  
  
  w0_prefix_sum <- c(0, cumsum(w0_sorted))  
  
  
  result <- 0
  i <- 1
  while (i <= n) {
    
    j <- i
    while (j <= n && s_sorted[j] == s_sorted[i]) j <- j + 1
    j <- j - 1  
    
    # s[i] > s[j] part
    if (i > 1) {
      result <- result + sum(w1_sorted[i:j]) * w0_prefix_sum[i]
    }
    
    # s[i] = s[j] part
    equal_w1_sum <- sum(w1_sorted[i:j])
    equal_w0_sum <- sum(w0_sorted[i:j])
    result <- result + 0.5 * equal_w1_sum * equal_w0_sum
    
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
  
  
  sum_auc=(sum_auc )/(sum(pi1)*sum(pi0)) #/P(Y=1)*P(Y=0) 
  return(sum_auc)
}


weighted_sum_kernel <- function(hb,s, w1, w0) {
  n <- length(s)
  
  diff_mat <- outer(s, s, "-") / hb
  
  # Phi matrix
  Phi_mat <- pnorm(diff_mat)
  
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
  sum_auc=(sum_auc )/(sum(pi1)*sum(pi0)) 
  return(sum_auc)
}

library(Rcpp)
sourceCpp("cpp/der.cpp")

sourceCpp("cpp/weighted_sum.cpp")
## IF kernel source part
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
  
  IF_1=scale( IF_1,center=T,scale=F)
  
  
  IF_2=IF_1-drop(est_f_no_sum(data = data,alpha =alpha ,nu=nu)%*%t(solve(A_total))%*%drop(der_AUC_alpha))
  
  return((IF_2))
  #return(mean(IF_2^2))
}
# derivative
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



## IF kernel source
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
  
  IF_1=scale( IF_1,center=T,scale=F)
  
  IF_2=IF_1-drop(est_f_no_sum(data = data,alpha =alpha ,nu=nu)%*%t(solve(A_total))%*%drop(der_AUC_alpha))#part 3_1
  
  return((IF_2))
  
}
# derivative
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



#----------------------auc-boot----------------

light_auc<- function(y, score) {
  r <- rank(score)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
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



