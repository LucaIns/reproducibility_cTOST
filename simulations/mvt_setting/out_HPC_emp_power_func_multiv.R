rm(list=ls())
# To print the command and its result in the output.
options(echo=TRUE)
# simulation inputs from bash
args <- commandArgs(trailingOnly = TRUE)
print(args)
sim_setting = as.numeric(args[1]) 
B = as.numeric(args[2]) # 10000 # 
load_dir       = paste0("./sim_HPC_out/sim_setting_", sim_setting, "/tmp")
destinationDir = paste0("./sim_HPC_out/sim_setting_", sim_setting)
##############################################
alpha = 0.05
delta = cte = log(1.25)
theta_= seq(0, 1.2, length=30)
if ((sim_setting==1) || (sim_setting==3)) {
  nu_ = 30    
} else if ((sim_setting==2) || (sim_setting==4)) {
  nu_ = 40
}
sigma1_ = c(0.08, 0.12, 0.16, 0.08, 0.08, 0.12)
sigma2_ = c(0.08, 0.12, 0.16, 0.12, 0.16, 0.16)
variab_=c("small","medium","large","mix_small","mix_medium","mix_large")
rho_ = c(0, 0.5, 0.9)
####################################
pow_emp_TOST_B_ =
  pow_emp_alphaTOST_B_ =
  pow_emp_xTOST_B_ = 
  array(NA,dim=c(length(nu_),length(variab_),length(rho_),length(theta_),B),
        dimnames=list(paste0("nu=",nu_),
                      paste0("variab=",variab_),
                      paste0("rho=",rho_),
                      paste0("theta=",theta_),
                      paste0("b=",1:B)))
ids_na = NULL
for(i in 1:B){
  if ((i%%1000) == 0){
    print(i)    
  }
  result = tryCatch({
      load(paste0(load_dir,"/",i))
  }, error = function(error_condition) {
      err=1
  }, warning = function(error_condition) {
      err=1})
  if(any(result==1)){
    ids_na = c(ids_na, i)
    next
  }
  if (any(is.na(pow_emp_alphaTOST_))){
    ids_na = c(ids_na, i)
  }
  pow_emp_TOST_B_[,,,,i] = pow_emp_TOST_
  pow_emp_alphaTOST_B_[,,,,i] = pow_emp_alphaTOST_
  pow_emp_xTOST_B_[,,,,i] = pow_emp_xTOST_
}
print(ids_na)
pow_emp_TOST = apply(pow_emp_TOST_B_,c(1,2,3,4),mean,na.rm=T)
pow_emp_alphaTOST = apply(pow_emp_alphaTOST_B_,c(1,2,3,4),mean,na.rm=T)
pow_emp_xTOST = apply(pow_emp_xTOST_B_,c(1,2,3,4),mean,na.rm=T)
save(pow_emp_TOST, file=paste0(destinationDir, "/pow_emp_TOST.rdata"))
save(pow_emp_alphaTOST, file=paste0(destinationDir, "/pow_emp_alphaTOST.rdata"))
save(pow_emp_xTOST, file=paste0(destinationDir, "/pow_emp_xTOST.rdata"))
##### end out HPC multiv ####