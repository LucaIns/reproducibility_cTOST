rm(list=ls())
## (Hyper)-parameters
# level
alpha = 0.05
# BE cte
delta = log(1.25)
# dof
nu_ = c(20,40,60,80)
# diff. BE parameter
theta_=c(0,delta)
# std. error
sigma_ = seq(0.01,0.2,l=1000)
# MC number of iterations
B = 1e5

## Recombining
# Directory
load_dir = "sim_HPC_out_curves_sigma/"
# To store results before computing mean
BE_TOST_B_ =
  BE_xTOST_no_B_ =
  BE_xTOST_off_B_ =
  BE_aTOST_B_ =
  array(NA,dim=c(length(nu_),length(sigma_),
                 length(theta_),B),
        dimnames=list(paste0("nu=",nu_),
                      paste0("sigma=",sigma_),
                      paste0("theta=",round(theta_,4)),
                      paste0("b=",1:B)))
incomplete_i=c()
# Recombining
for(i in 1:B){
  #
  if ((i%%100) == 0){
    print(i)
  }
  #
  res=tryCatch({load(paste0(load_dir,i))},error=function(e){e})
  if(inherits(res,"error")){
    incomplete_i = c(incomplete_i,i)
    next
  }  #
  BE_TOST_B_[,,,i] = BE_TOST_
  BE_xTOST_no_B_[,,,i] = BE_xTOST_no_
  BE_xTOST_off_B_[,,,i] = BE_xTOST_off_
  BE_aTOST_B_[,,,i] = BE_aTOST_
}

pow_TOST_curves_sigma = apply(BE_TOST_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSTno_curves_sigma = apply(BE_xTOST_no_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSToff_curves_sigma = apply(BE_xTOST_off_B_,c(1,2,3),mean,na.rm=T)
pow_aTOST_curves_sigma = apply(BE_aTOST_B_,c(1,2,3),mean,na.rm=T)

save(pow_TOST_curves_sigma, file="pow_TOST_curves_sigma.rdata")
save(pow_xTOSTno_curves_sigma, file="pow_xTOSTno_curves_sigma.rdata")
save(pow_xTOSToff_curves_sigma, file="pow_xTOSToff_curves_sigma.rdata")
save(pow_aTOST_curves_sigma, file="pow_aTOST_curves_sigma.rdata")



