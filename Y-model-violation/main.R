
#-----------------prepare-----------------------------------------


source(file = "main-fun-parallel.R")
source(file="parallel_other_fun.R")
sourceCpp(file = "cpp/der.cpp")
sourceCpp(file = "cpp/weighted_sum.cpp")
source(file="tune-param.R")

#------------------true target-------------------

# main effect of W on Y
theta_t <- c(
  0,          
  -log(1.8),
  log(2.5),
  log(1.4),
  log(1.8)
)
# Interaction effects between Z and selected X variables
gamma_xz <- c(0.2, -0.2, 0.15)

# Target outcome model parameters:
# includes main effects of W and X-Z interaction terms

yparameter_t <- list(
  beta = theta_t,
  gamma_xz = gamma_xz,
  include_xz_interaction = TRUE
)
zparam = list(beta= c(0.4,0.5,-0.2,0.3), sigma_z=1)


# calibrate target intercept
theta_t0<-calibrate_T_intercept(5000000,yparameter_t = yparameter_t,prevalence = 0.1,zparameter = zparam)

yparameter_t$beta[1]<-theta_t0


# true value
True_T_param<-true_T_est(nt = 5000000,yparameter_t = yparameter_t,zparameter = zparam)



#----------------------source-true(OM2)--------------------
yparameter_s= yparameter_t
True_S_param<-true_S_est(nt=2000,ns=5000000,yparameter_s =yparameter_s ,
           xi=c(0.5,-0.5,0.3),
           zparameter=zparam)
True_T_param
True_S_param
#-------------------模拟-----------

nt=1000
ns=500

reslist=list()



# OM1

t1=Sys.time()

res_single<-main_simu2_parallel(nco=10,N=1000,nt=nt,ns=ns,yparameter_s = yparameter_s,
                                yparameter_t=yparameter_t,
                                xi=c(0,0,0),
                                zparameter = zparam)
t2=Sys.time()
t2-t1
reslist[[1]]=res_single
gc()

# OM2

t1=Sys.time()

res_single<-main_simu2_parallel(nco=10,N=1000,nt=nt,ns=ns,yparameter_s = yparameter_s,
                                yparameter_t=yparameter_t,xi=c(0.5,-0.5,0.3),
                                zparameter = zparam)
t2=Sys.time()
t2-t1
reslist[[2]]=res_single
gc()


save(reslist,zparam,evalua,yparameter_t,yparameter_s,True_T_param,True_S_param,file="res/summary_navie_ci_ns500_OM.rdata")



