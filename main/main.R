
#-----------------source code-----------------------------------------


source(file = "main-fun-parallel.R")
source(file="parallel_other_fun.R")
sourceCpp(file = "cpp/der.cpp")
sourceCpp(file = "cpp/weighted_sum.cpp")
#Sys.setenv(OMP_NUM_THREADS=4)


#------------------true-------------------

betat=c(-3.05,-log(1.8),log(2.5),log(1.4),log(1.8))
zparam = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=1)

## auc true value (auc for full model and base model; auc difference)
auc_true_xz<-function(nt,betat,mu_xt=rep(0,3),sigma_xt=0.5^ abs(outer(1:3, 1:3, "-")), 
                      zparameter = zparam,
                      threshold_p=0) {
  set.seed(2)
  p1 <- length(mu_xt)
  
  # -------------------
  # Target: generate (X_{-p}, latent X_p) together
  # -------------------
  # continuous part
  xt_cont <- mvrnorm(nt, mu = mu_xt, Sigma = sigma_xt)  # last column is latent X_p
  # discrete
  xp_t <- as.numeric(xt_cont[,p1] > threshold_p)
  
  xt <- cbind(xt_cont[,-p1], xp_t)
  
  # Z | X
  beta_z <- zparameter$beta
  sigma_z <- zparameter$sigma_z
  
  zt <- cbind(1, xt) %*% beta_z[1:(p1+1)] + rnorm(nt,0,sigma_z)
  
  
  # Y | X, Z
  prob_yt <- 1/(1+exp(-cbind(1, xt, zt) %*% betat))
  yt <- rbinom(nt,1,prob_yt)
  
  
  yt_model0 <- glm(yt ~ xt, family = binomial)
  roc_obj0=roc(yt,yt_model0$fitted.values)
  print(mean(yt))
  
  yt_model1 <- glm(yt ~ xt+zt, family = binomial)
  roc_obj1=roc(yt,yt_model1$fitted.values)
  
  auc0=as.numeric(auc(roc_obj0))
  auc1=as.numeric(auc(roc_obj1))
  aucdiff=auc1-auc0
  return(c(auc1,auc0,aucdiff) )
  
  
}
t1=Sys.time()
auc_true_est=auc_true_xz(nt = 5000000,betat = betat)
t2=Sys.time()
t2-t1



#-------------------simulation-----------

nt=1000
ns=500

betat=c(-3.05,-log(1.8),log(2.5),log(1.4),log(1.8))
zparam = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=1)

reslist=list()


rhox=c(0.1, -0.2, 0.3, 0.2)

rho_list_diffz=list(rep(0,5),rep(0,5),
                    c(rhox,  0),
                    c(rhox,  -0.3),c(rhox,   0.3))


# C1
rho_beta=rho_list_diffz[[1]]
betas=betat*(1+rho_beta)
t1=Sys.time()
# nco: core number
res_single<-main_simu2_parallel(nco=9,N=1000,nt=nt,ns=ns,betat=betat,betas=betas,xi=c(0,0,0),
                                zparameter = zparam)
t2=Sys.time()
t2-t1
reslist[[1]]=res_single
gc()

# C2
rho_beta=rho_list_diffz[[2]]
betas=betat*(1+rho_beta)
t1=Sys.time()

res_single<-main_simu2_parallel(nco=9,N=1000,nt=nt,ns=ns,betat=betat,betas=betas,xi=c(0.5,-0.5,0.3),
                                zparameter = zparam)
t2=Sys.time()
t2-t1
reslist[[2]]=res_single
gc()

# C3 C4a C4b


for(i in 3:length(rho_list_diffz) ){
  rho_beta=rho_list_diffz[[i]]
  betas=betat*(1+rho_beta)
  t1=Sys.time()
  
  res_single<-main_simu2_parallel(nco=9,N=1000,nt=nt,ns=ns,betat=betat,betas=betas,xi=c(0.5,-0.5,0.3),
                                  zparameter = zparam)
  t2=Sys.time()
  cat("\n",t2-t1)
  reslist[[i]]=res_single
  gc()
}

save(reslist,rho_list_diffz,zparam,evalua,betat,auc_true_est,file="res/summary_navie_ci_ns500.rdata")
