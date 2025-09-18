rm(list=ls())
## (Hyper)-parameters
# level
alpha = 0.05
# BE cte
delta = log(1.25)
# diff. BE parameter
gamma_=c(1,1.25)
theta_=log(gamma_)
# dof
nu_ = c(seq(5,80,by=1),100,250,500,1000)
# std. error
sigma_ = seq(0.005,0.3,by=0.005)

## MC params
B=10^5

## Recombining
# Directory
load_dir = "out_tiles/"
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
#
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
  }
  #
  BE_TOST_B_[,,,i] = BE_TOST_
  BE_xTOST_no_B_[,,,i] = BE_xTOST_no_
  BE_xTOST_off_B_[,,,i] = BE_xTOST_off_
  BE_aTOST_B_[,,,i] = BE_aTOST_
}
#
pow_TOST_tiles = apply(BE_TOST_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSTno_tiles = apply(BE_xTOST_no_B_,c(1,2,3),mean,na.rm=T)
pow_xTOSToff_tiles = apply(BE_xTOST_off_B_,c(1,2,3),mean,na.rm=T)
pow_aTOST_tiles = apply(BE_aTOST_B_,c(1,2,3),mean,na.rm=T)
#
save(pow_TOST_tiles, file="pow_TOST_tiles.rdata")
save(pow_xTOSTno_tiles, file="pow_xTOSTno_tiles.rdata")
save(pow_xTOSToff_tiles, file="pow_xTOSToff_tiles.rdata")
save(pow_aTOST_tiles, file="pow_aTOST_tiles.rdata")







