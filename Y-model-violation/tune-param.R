#-------------------------test function--------------------------


# target intercept
calibrate_T_intercept <- function(nt,yparameter_t,prevalence,
                                  mu_xt=rep(0,3), sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), 
                                  zparameter = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=0.7, include_interaction=FALSE),
                                  threshold_p=0,lower = -20, upper = 20) {
  set.seed(1)
  p1 <- length(mu_xt)
  expit <- function(u) 1 / (1 + exp(-u))
  
  # -------------------
  # Helper: linear predictor for Y | X, Z
  # -------------------
  yparameter_t$beta[1]<-0
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
  sigma_z <- zparameter$sigma
  if(zparameter$include_interaction){
    #
  } else {
    zt <- cbind(1, xt) %*% beta_z[1:(p1+1)] + rnorm(nt,0,sigma_z)
  }
  
  
  # -------------------
  # Target: Y | X, Z
  # -------------------
  
  eta_yt <- eta_y_fun(xt, zt, yparameter_t)
  
  uniroot(function(intercept0){mean(expit(intercept0+eta_yt))-prevalence}, lower = lower, upper = upper)$root
  
  
}


# target betaT(pseudo parameter) true value


true_T_est <- function(nt,yparameter_t,mu_xt=rep(0,3),
                       sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")),
                       zparameter,threshold_p=0) {
  set.seed(1)
  p1 <- length(mu_xt)
  expit <- function(u) 1 / (1 + exp(-u))
  
  # -------------------
  #  linear predictor for Y | X, Z
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
  # Working models under target population
  # -------------------
  yt_model0 <- glm(yt ~ xt, family = binomial)
  psi_t_dagger <- coef(yt_model0)
  
  yt_model1 <- glm(yt ~ xt + zt, family = binomial)
  beta_t_dagger <- coef(yt_model1)
  
  
  # -------------------
  # AUC truth based on working scores
  # -------------------
  score0 <- drop(predict(yt_model0, type = "link"))
  score1 <- drop(predict(yt_model1, type = "link"))
  
  roc_obj0 <- pROC::roc(yt, score0, quiet = TRUE)
  roc_obj1 <- pROC::roc(yt, score1, quiet = TRUE)
  
  auc0 <- as.numeric(pROC::auc(roc_obj0))
  auc1 <- as.numeric(pROC::auc(roc_obj1))
  aucdiff <- auc1 - auc0
  
  
  
  return(list(
    auc = c(auc1 = auc1, auc0 = auc0, aucdiff = aucdiff),
    beta_t_dagger = beta_t_dagger,
    psi_t_dagger = psi_t_dagger,
    prevalence = mean(yt)
  ))
}

true_S_est <- function(nt,ns,yparameter_s,mu_xt=rep(0,3),
                       sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")),xi,
                       zparameter = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=0.7, include_interaction=FALSE),
                       threshold_p=0) {
  set.seed(1)
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
  sigma_z <- zparameter$sigma
  
  
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
  
  # -------------------
  # Working models under source population
  # -------------------
  
  
  ys_model1 <- glm(ys ~ xs + zs, family = binomial)
  beta_s_dagger <- coef(ys_model1)
  
  
  
  
  
  return(list(
    beta_s_dagger = beta_s_dagger,
    prevalence = mean(ys)
  ))
}
