library(pracma)
library(sandwich)
library(nloptr)
library(pROC)
library(withr)
library(MASS)



data.gen2 <- function(ns, nt, 
                      yparameter_s,yparameter_t,
                      mu_xt=rep(0,3), sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), xi, 
                      zparameter,threshold_p=0) {
  
  p1 <- length(mu_xt)
  expit <- function(u) 1 / (1 + exp(-u))
  
  # -------------------
  # Helper: linear predictor for Y | X, Z
  # -------------------
  eta_y_fun <- function(x, z, yparameter) {
    beta <- yparameter$beta
    eta <- drop(cbind(1, x, z) %*% beta)
    
    if (isTRUE(yparameter$include_xz_interaction)) {
      gamma_xz <- yparameter$gamma_xz
      eta <- eta + drop((x * drop(z)) %*% gamma_xz)
    }
    
    return(eta)
  }
  # -------------------
  # Target: generate (X_{-p}, latent X_p) together
  # -------------------
  xt_cont <- mvrnorm(nt, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  
  xp_t <- as.numeric(xt_cont[,p1] > threshold_p)
  xt <- cbind(xt_cont[,-p1], xp_t)
  
  # Z | X
  beta_z <- zparameter$beta
  sigma_z <- zparameter$sigma_z
  zt <- cbind(1, xt) %*% beta_z[1:(p1+1)] + rnorm(nt,0,sigma_z)
  
  

  # -------------------
  # Target: Y | X, Z
  # -------------------
  eta_yt <- eta_y_fun(xt, zt, yparameter_t)
  prob_yt <- expit(eta_yt)
  yt <- rbinom(nt, 1, prob_yt)

  # -------------------
  # Source: tilting for X
  # -------------------
  
  # Source covariates are generated under the density ratio model
  # f_S(X) \propto f_T(X) exp(-X^T xi). Since the continuous latent
  # variables are Gaussian, the tilting for the continuous components
  # can be incorporated through a mean shift in the proposal distribution.
  # The remaining tilting for the thresholded binary covariate is handled
  # by rejection sampling.
  
  # Mean shift induced by exponential tilting of the continuous components.
  # The last latent variable is also shifted because it is correlated with
  # the continuous covariates.
  adjust_p=drop(sigma_xt[p1,-p1]) %*% xi[-p1]
  mu_proposal <- mu_xt - c(drop(sigma_xt[-p1,-p1] %*% xi[-p1]),adjust_p)
  
  # Rejection sampling accounts for the remaining tilting term
  # exp(-xi_p X_p) for the binary covariate.
  # Compute the acceptance probability for the binary covariate component.
  M = max(1, exp(-xi[p1]))
  p_accept=( mean(xp_t)*exp(-xi[p1])+1-mean(xp_t) )/M
  n_temp=ceiling(1.5 * ns / p_accept)
  
  # Generate candidate samples from the proposal distribution.
  x_cont_candidate  <- mvrnorm(n_temp, mu = mu_proposal, Sigma = sigma_xt)  # includes latent 
  # Dichotomize the latent variable to obtain candidate binary covariates.
  xp_candidate <- as.numeric(x_cont_candidate[,p1] >  threshold_p)
  
  # Rejection sampling for the binary component under exponential tilting.
  
  u <- runif(n_temp)
  accept <- u < exp(-xi[p1]*xp_candidate)/M
  x_cont_s <- x_cont_candidate[accept,,drop=FALSE]
  xp_s <- xp_candidate[accept]
  
  
  xs <- cbind(x_cont_s[,-p1,drop=FALSE], xp_s)
  xs=xs[1:ns,]
  # Z | X 
  zs <- cbind(1, xs) %*% beta_z[1:(p1+1)] + rnorm(ns,0,sigma_z)
  
  # -------------------
  # Source: Y | X, Z
  # -------------------
  eta_ys <- eta_y_fun(xs, zs, yparameter_s)
  prob_ys <- expit(eta_ys)
  ys <- rbinom(ns, 1, prob_ys)
  #-----------------summary-----------
  yt_xt_model <- glm(yt ~ xt, family = binomial)
  
  
  pred_p=drop(cbind(1,xt)%*%coef(yt_xt_model))
  pred_p=1/(1+exp(-pred_p))
  
  A_ma=t(cbind(1,xt))%*%(cbind(1,xt)*pred_p*(1-pred_p))/nrow(xt)
  
  IF_psi=(cbind(1,xt)*(yt-pred_p)) %*% t(solve(A_ma))
  var_nu=cov(cbind(1,xt,IF_psi))/nrow(xt)
  

  #-----------results-------------
  return(list(xs=xs,zs=zs,ys=ys,mut=colMeans(cbind(1,xt)),
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



