evalua<-function(res,p=4,q=1,beta_t){
  
  # targetbeta=colMeans(result$beta_cmle)
  
  # tauc[1]=mean(result$auchat_cmle[,1])
  # tauc[2]=mean(result$auchat_reduce[,1])
  # tauc[3]=mean(result$auchat_dif[,1])
  ind_t=c((2*p+q+1):(3*p+q),p+q)
  total_comp<-data.frame(c("S1","","","",""),
                         c("$\\beta_{T,X_0}$","$\\beta_{T,X_1}$","$\\beta_{T,X_2}$","$\\beta_{T,X_3}$","$\\beta_{T,Z}$"),              
                         beta_t,
                         
                         colMeans(res$res_trans[,ind_t])-beta_t,
                         sqrt(diag(var(res$res_trans[,ind_t]))),
                         colMeans(res$sd_trans[,ind_t]),
                         colMeans(t( (t(res$res_trans[,ind_t])-beta_t)^2)),
                         rowMeans(abs(t(res$res_trans[,ind_t])-beta_t)/t(res$sd_trans[,ind_t])<qnorm(0.975)),
                         colMeans(res$res_mle)-beta_t,
                         sqrt(diag(var(res$res_mle))),
                         colMeans(res$sd_mle),
                         colMeans(t( (t(res$res_mle)-beta_t)^2)),
                         rowMeans(abs(t(res$res_mle)-beta_t)/t(res$sd_mle)<qnorm(0.975))
                         #c(rowMeans(abs(t(result$beta_cmle[,1:5])-targetbeta[1:5])/t(result$asyd_cml)<qnorm(0.975)),0,mean(abs(result$auchat_cmle[,1]-tauc[1])/sqrt(result$auchat_cmle[,2])<qnorm(0.975) ),mean(abs(result$auchat_dif[,1]-mean(result$auchat_dif[,1]))/sqrt(result$auchat_dif[,2])<qnorm(0.975) ) )
                         
  )
  auctrue=True_T_param$auc
  auc_xz_true_est=auctrue[1]
  auc_x_true_est=auctrue[2]
  auc_diff_true_est=auctrue[3]
  other_metric<-data.frame(c("","",""),
                           
                           c("$\\AUC$  w/ $Z$","$\\AUC$  w/o $Z$","$\\Delta{\\AUC}$"),
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
