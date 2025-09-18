rm(list=ls())
## (Hyper)-parameters
# level
alpha = 0.05
# BE cte
delta = log(1.25)
# diff. BE parameter
theta_=seq(0,delta*1.2,by=delta/28)
# dof
nu_ = c(5,10,15,20,25,30,40,50,60)
# std. error
sigma_ = c(8, 12, 16)/100
## MC params
B=10^5

## Recombining
# Directory
load_dir = "sim_HPC_out_power_curves/"
# To store results before computing mean
BE_TOST_B_ =
  BE_xTOST_no_B_ =
  BE_xTOST_off_B_ =
  BE_xTOST_boot_B_ =
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
  BE_xTOST_boot_B_[,,,i] = BE_xTOST_boot_
  BE_aTOST_B_[,,,i] = BE_aTOST_
}
#
pow_TOST_power_curves = apply(BE_TOST_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSTno_power_curves = apply(BE_xTOST_no_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSToff_power_curves = apply(BE_xTOST_off_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSTboot_power_curves = apply(BE_xTOST_boot_B_,c(1,2,3),mean,na.rm=T)
pow_aTOST_power_curves = apply(BE_aTOST_B_,c(1,2,3),mean,na.rm=T)
#
save(pow_TOST_power_curves, file="pow_TOST_power_curves.rdata")
save(pow_xTOSTno_power_curves, file="pow_xTOSTno_power_curves.rdata")
save(pow_xTOSToff_power_curves, file="pow_xTOSToff_power_curves.rdata")
save(pow_xTOSTboot_power_curves, file="pow_xTOSTboot_power_curves.rdata")
save(pow_aTOST_power_curves, file="pow_aTOST_power_curves.rdata")



