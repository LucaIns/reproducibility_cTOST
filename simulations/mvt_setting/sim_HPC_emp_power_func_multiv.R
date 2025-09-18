rm(list=ls())

# source the functions to use
source("aux_fun/atost_mvt.R")
source("aux_fun/ctost_univ.R")
source("aux_fun/ctost_mvt.R")

# print the command and its result in the output.
options(echo=TRUE)
clusterUse=T
# clusterUse=F
library(tictoc)
##########################################################
# number of simulations and arrays
if (clusterUse==T) {
  # simulation inputs from bash
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  sim_setting = as.numeric(args[1])
  B = n_simu = as.numeric(args[2])
  n_array = as.numeric(args[3])
  # create matrix of indices
  ind_mat <- matrix(1:n_simu, nr = n_array, byr = T)
  print(ind_mat)
  # get slurm array id and convert to numeric
  id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")
  id_slurm <- as.numeric(id_slurm)
  print("id_slurm: ")
  print(id_slurm)
  # define id of simu to be run on array
  Rw <- ind_mat[id_slurm, ]
  # create the working directory
  destinationDir <- paste0("sim_HPC_out/sim_setting_", sim_setting, "/tmp") 
  if(!dir.exists(destinationDir)) dir.create(destinationDir,recursive=T)
} else if (clusterUse==F) {
  destinationDir=""
  B  = 100
  Rw = 1:B
  sim_setting = 1
}
print(Rw)
############################################################
alpha = 0.05
delta = cte = log(1.25)
theta_= seq(0, 1.2, length=30)
# number of MC for power function (mvt alpha-TOST)
B_aTOST = 5*10^4
# cTOST options
max_iter = 10
tolpower=0.0001
# sim parameters
sigma1_ = c(0.08, 0.12, 0.16, 0.08, 0.08, 0.12)
sigma2_ = c(0.08, 0.12, 0.16, 0.12, 0.16, 0.16)
variab_=c("small","medium","large","mix_small","mix_medium","mix_large")
rho_ = c(0, 0.5, 0.9)
if ((sim_setting==1) || (sim_setting==3)) {
  nu_ = 30    
} else if ((sim_setting==2) || (sim_setting==4)) {
  nu_ = 40
}
# initialization
tt_nu_=rep(NA,length(nu_))
tt_variab_=rep(NA,length(variab_))
tt_rho_=rep(NA,length(rho_))
tt_Rw_=rep(NA,length(Rw))
tic()
cat("Start -", date(),"\n\n\n")
for(rw in Rw){
  if(!any(dir(destinationDir)==rw)){ # if this solution is not already present
    tic()
    cat("start Rw:",which(Rw==rw),"/",length(Rw),"\n")
    pow_emp_TOST_ =
      pow_emp_alphaTOST_ =
      pow_emp_xTOST_ = array(NA,dim=c(length(nu_),length(variab_),
                                      length(rho_),length(theta_)),
                             dimnames=list(paste0("nu=",nu_),
                                           paste0("variab=",variab_),
                                           paste0("rho=",rho_),
                                           paste0("theta=",theta_)))
    for(u in 1:length(nu_)){
      tic()
      if(u>1){
        cat("Rw:",which(Rw==rw),"/",length(Rw),"\n")
        cat("\tstart nu_:",u,"/",length(nu_),"\n")
      }
      counter_pb = 0
      pb = txtProgressBar(min = 0, max = length(variab_)*length(rho_), style = 3)
      nu = nu_[u]
      for(v in 1:length(variab_)){
        if ((sim_setting==1) || (sim_setting==2)) {
          sigmas = c(sigma1_[v], sigma2_[v])                         # bivariate
        } else if ((sim_setting==3) || (sim_setting==4)) {
          sigmas = c(sigma1_[v], sigma1_[v], sigma2_[v], sigma2_[v]) # 4-variate
        }
        for(r in 1:length(rho_)){ 
          rho = rho_[r]
          if ((sim_setting==1) || (sim_setting==2)) {
            Sigma = build_Sigma(sigmas=sigmas,rho)
          } else if ((sim_setting==3) || (sim_setting==4)) {
            Sigma = build_Sigma(sigmas=sigmas,rho,type="ar1")
          }
          set.seed(rw)
          Sigma.hat = sim_Sigma.hat_mv(Sigma,nu)
          # alpha-TOST
          out_alphaTOST_sol = get_alpha_TOST_MC_mv(
            alpha=alpha,
            Sigma=Sigma.hat,
            nu=nu,
            delta=cte,
            B=B_aTOST)
          out_alphaTOST = out_alphaTOST_sol$min
          # cTOST
          out_xTOST_sol = get_ctost_mvt(alpha, Sigma.hat, rep(cte, length(sigmas)), 
            max_iter=max_iter,
            theta=NULL, 
            tolpower=tolpower)
          out_xTOST = out_xTOST_sol$c_of_0
          # get argsup at the population level
          # alpha-TOST
          out_alphaTOST_sol_pop = get_alpha_TOST_MC_mv(
                      alpha=alpha,
                      Sigma=Sigma,
                      nu=nu,
                      delta=cte,
                      B=B_aTOST)
          out_alphaTOST_pop = out_alphaTOST_sol_pop$min
          theta_sup_alphaTOST = find_sup_x(out_alphaTOST_pop, Sigma, delta)
          theta_sup_TOST = find_sup_x(alpha, Sigma, delta)
          # cTOST
          out_xTOST_sol_pop = get_ctost_mvt(alpha, Sigma, rep(cte, length(sigmas)),
            max_iter=max_iter,
            theta=NULL, 
            tolpower=tolpower)
          # TOST
          theta_sup_xTOST = find_sup_ctost(Sigma, rep(cte, length(sigmas)), out_xTOST_sol_pop$c_of_0)
          for(k in 1:length(theta_)){
            # method-specific theta seq
            # TOST
            theta_TOST = theta_sup_TOST * theta_[k]
            set.seed(rw+B)
            theta.hat_TOST = sim_theta.hat_mv(theta_TOST,Sigma)
            # alpha-TOST
            theta_aTOST = theta_sup_alphaTOST * theta_[k] 
            set.seed(rw+B)
            theta.hat_aTOST = sim_theta.hat_mv(theta_aTOST,Sigma)
            # cTOST
            theta_xTOST = theta_sup_xTOST * theta_[k] 
            set.seed(rw+B)
            theta.hat_xTOST = sim_theta.hat_mv(theta_xTOST,Sigma)
            # construct CI
            # TOST
            ci_TOST = ci_k(theta.hat_TOST, Sigma.hat, nu, alpha)
            pow_BE_TOST = BE_mv(ci_TOST,cte)
            pow_emp_TOST_[u,v,r,k] = pow_BE_TOST
            # alpha-TOST
            ci_alphaTOST = ci_k(theta.hat_aTOST, Sigma.hat, nu, out_alphaTOST)
            pow_BE_alphaTOST = BE_mv(ci_alphaTOST,cte)
            pow_emp_alphaTOST_[u,v,r,k] = pow_BE_alphaTOST
            # cTOST
            pow_BE_xTOST = prod(abs(theta.hat_xTOST) < out_xTOST)
            pow_emp_xTOST_[u,v,r,k] = pow_BE_xTOST
          } # end theta
          counter_pb = counter_pb + 1
          setTxtProgressBar(pb, counter_pb)
        } # end rho
      } # end sigma
      tt_nu=toc(quiet=T)
      tt_nu_[u]=tt_nu$toc-tt_nu$tic
      if(u==length(nu_)){
        cat("\n\tend nu_:",u,"/",length(nu_),"- ",tt_nu$callback_msg,"\n")
      }else{
        cat("\n\tend nu_:",u,"/",length(nu_),"- ",tt_nu$callback_msg,"\n\n")
      }
    } # end nu
    tt_Rw=toc(quiet=T)
    tt_Rw_[which(Rw==rw)]=tt_Rw$toc-tt_Rw$tic
    cat("end Rw:",which(Rw==rw),"/",length(Rw),"- ",tt_Rw$callback_msg,"\n\n")
    # save param
    if(clusterUse==T){
      save(
        pow_emp_TOST_,
        pow_emp_alphaTOST_,
        pow_emp_xTOST_,
        file=paste0(destinationDir,"/",rw))
    } 
  } # end if
} # end rw
tt = toc(quiet=T)
t = tt$toc-tt$tic
cat("End -", tt$callback_msg, "-", date(),"\n")
##### end sim HPC multiv ####