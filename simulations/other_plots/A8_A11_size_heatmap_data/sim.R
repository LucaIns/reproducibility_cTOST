rm(list=ls())
options(scipen = 999)
clusterUse=T
load("xtost_offline_grid.rda")
#### Functions & libraries ####
library(PowerTOST)
library(tictoc)
alphahat.fun = function(sigma.hat,alpha,delta,nu,tol=1e-8){
  K = 10000
  alpha.k = c(alpha,rep(NA,K-1))
  for(k in 2:K){
    tval       = qt(1 - alpha.k[k-1], df = nu)
    delta1     = (2*delta)/sigma.hat
    delta2     = 0
    R          = (delta*sqrt(nu))/(tval*sigma.hat)
    omega      = PowerTOST::OwensQ(nu, -tval, delta2, 0, R)-
      PowerTOST::OwensQ(nu,  tval, delta1, 0, R)
    alpha.k[k] = min(c(alpha + alpha.k[k-1] - omega,0.5))
    if(abs(alpha.k[k]-alpha.k[k-1])<tol){break}
  }
  # out
  ifelse(k==K,NA,alpha.k[k])
}
# Function to simulate data for boot xTOST
simulate_data = function(mu, sigma, nu, seed = 18){
  set.seed(seed)
  theta_hat = rnorm(1, mean = mu, sd = sigma)
  sig_hat = sqrt(rchisq(1, df = nu)/nu*sigma^2)
  list(theta_hat = theta_hat, sig_hat = sig_hat)
}
# Function to compute the margins for xTOST
get_c_of_0 = function(delta, sig_hat, alpha, B = 1000, tol = 10^(-8)){
  c0 = delta
  alpha0 = alpha
  
  # Sequence of c's
  m = 30
  cte_vect = rep(NA, m)
  
  # Initial approx
  c_init = c0 - sig_hat*qnorm(1 - alpha0)
  if (c_init < 0){
    c_init = c0/2
  }
  cte_vect[1] = c_init
  
  # Start newton raphson
  for (i in 1:(B-1)){
    delta = (pnorm((c0 + cte_vect[i])/sig_hat) - pnorm((c0 - cte_vect[i])/sig_hat) - alpha0)*sig_hat /
      (dnorm((c0 + cte_vect[i])/sig_hat) + dnorm((c0 - cte_vect[i])/sig_hat))
    cte_vect[i+1] = cte_vect[i] - delta
    
    if (cte_vect[i+1] < tol || cte_vect[i+1] > c0){
      cte_vect[i+1] = c0 - tol
    }
    
    if (abs(delta) < tol){
      out = list(c = cte_vect[i+1], converged = TRUE, iter = i+1)
      break
    }
  }
  
  if(i < B-1){
    size_NR = size_xTOST(sig_hat=sig_hat, delta=c0, delta_star=out$c)
    out = list(c = out$c, size=size_NR, converged = out$converged, iter = out$iter)
  }else{
    size_NR = size_xTOST(sig_hat=sig_hat, delta=c0, delta_star=cte_vect[B])
    out = list(c = cte_vect[B], size=size_NR, converged = F, iter = B)
  }
  
  out
}
# Function to compute xTOST solutions
xtost = function(sig_hat, nu, alpha, delta, correction = "no", B = 10^5, theta_hat=NULL, seed = 85){
  if (!(correction %in% c("no", "bootstrap", "offline"))){
    stop("Finite sample correction method not implemented.")
  }
  if (correction == "bootstrap"){
    res = rep(NA, B)
    for (i in 1:B){
      dat = simulate_data(mu = delta, sigma = sig_hat,
                          nu = nu, seed = seed + i)
      c_0_hat = get_c_of_0(delta = delta, sig_hat = dat$sig_hat, alpha = alpha)
      res[i] = abs(dat$theta_hat) < c_0_hat$c
    }
    correct_alpha = 2*alpha - mean(res)
  }
  if (correction == "offline"){
    if(nu<=100){
      data_to_correct = xtost_offline_grid
      index_sigma = which.min(abs(data_to_correct$sigma - sig_hat))
      index_nu = which.min(abs(data_to_correct$nu - nu))
      correct_alpha = 2*alpha - data_to_correct$tier[index_nu, index_sigma]
    }else{
      correct_alpha = alpha # i.e. no correction
    }
  }
  if (correction == "no"){
    correct_alpha = alpha # i.e. no correction
  }
  c_0_hat = get_c_of_0(delta = delta, sig_hat = sig_hat, alpha = correct_alpha)
  return(c_0_hat)
  decision = abs(theta_hat) < c_0_hat$c
  ci_half_length = delta - c_0_hat$c
  ci = theta_hat + c(-1, 1) * ci_half_length
  out = list(decision = decision, ci = ci, theta_hat = theta_hat,
             sig_hat = sig_hat, nu = nu, alpha = alpha,
             c0 = c_0_hat$c,
             correction = correction,
             correct_alpha = correct_alpha,
             delta = delta, method = "x-TOST")
  class(out) = "tost"
  out
}
# Function to compute size of xTOST
size_xTOST = function(sig_hat, delta, delta_star, ...){
  power_xTOST(theta = delta, sig_hat = sig_hat, delta = delta_star)
}
# Function to compute power of xTOST
power_xTOST = function(theta, sig_hat, delta){
  pnorm((theta + delta)/sig_hat) - pnorm((theta - delta)/sig_hat)
}
#### End Functions & libraries ####
#### Simulations ####
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
## MC & Cluster
destinationDir <- "out_tiles"
if(!dir.exists(destinationDir)) dir.create(destinationDir)
if(clusterUse) jobId <- commandArgs(trailingOnly = TRUE) else jobId=1
jobId <- as.numeric(jobId)
#
B = 10^5
Rw = seq(as.numeric(jobId),B,5000)
#
tt_nu_=rep(NA,length(nu_))
tt_sigma_=rep(NA,length(sigma_))
tt_Rw_=rep(NA,length(Rw))
#
## Start simulations
cat("Start -", date(),"\n\n\n")
for(rw in Rw){
  # if(sum(rw==toeval)==1){ # if this solution is not already present
  tic()
  cat("start Rw:",which(Rw==rw),"/",length(Rw),"\n")
  # To store the results for each rw
  BE_TOST_ =
    BE_xTOST_no_ =
    BE_xTOST_off_ =
    BE_aTOST_ =
    array(NA,dim=c(length(nu_),length(sigma_),
                   length(theta_)),
          dimnames=list(paste0("nu=",nu_),
                        paste0("sigma=",sigma_),
                        paste0("theta=",round(theta_,4))))
  for(u in 1:length(nu_)){
    tic()
    nu = nu_[u]
    for(s in 1:length(sigma_)){
      sigma = sigma_[s]
      # Drawing sigma_hat
      sig_hat = sqrt(rchisq(1, df = nu)/nu*sigma^2)
      # 1. xTOST no correction
      c0_no_star = xtost(sig_hat = sig_hat, nu = nu, alpha = alpha, delta = delta, correction = "no", B = 10^5)
      # 2. xTOST offline solution
      c0_off_star = xtost(sig_hat = sig_hat, nu = nu, alpha = alpha, delta = delta, correction = "offline", B = 10^5)
      # 3. aTOST solution
      alpha_star = alphahat.fun(sigma.hat = sig_hat, alpha = alpha, delta = delta, nu = nu)
      for(k in 1:length(theta_)){
        theta = theta_[k]        
        # Drawing theta_hat
        theta_hat = rnorm(1, mean = theta, sd = sigma)
        # 0. TOST BE
        BE_TOST_[u,s,k] = abs(theta_hat) < (delta - qt(1 - alpha, df = nu)* sig_hat)
        # 1. xTOST no BE
        BE_xTOST_no_[u,s,k] = abs(theta_hat) < c0_no_star$c
        # 2. xTOST offline BE
        BE_xTOST_off_[u,s,k] = abs(theta_hat) < c0_off_star$c
        # 3. aTOST BE
        BE_aTOST_[u,s,k] = abs(theta_hat) < (delta - qt(1 - alpha_star, df = nu)* sig_hat)
      } # end theta_
    } # end sigma_
    tt_nu=toc(quiet=T)
    tt_nu_[u]=tt_nu$toc-tt_nu$tic
  } # end nu_
  tt_Rw=toc(quiet=T)
  tt_Rw_[which(Rw==rw)]=tt_Rw$toc-tt_Rw$tic
  cat("end Rw:",which(Rw==rw),"/",length(Rw),"- ",tt_Rw$callback_msg,"\n\n")
  # save param
  if(clusterUse==T){
    save(
      BE_TOST_,
      BE_xTOST_no_,
      BE_xTOST_off_,
      BE_aTOST_,
      file=paste0(destinationDir,"/",rw))
  }
} # end rw
mean(tt_Rw_)
cat("End -", date(),"\n")
#### End Simulations ####



