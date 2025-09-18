library(mvtnorm)
library(PowerTOST)
library(MCMCpack)
library(utils)

ar1_cor <- function(p, rho) {
  # construct correlation matrix
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) -
                    (1:p - 1))
  rho^exponent
}

ar1_cov <- function(sigmas, rho){
  p = length(sigmas)
  R = ar1_cor(p, rho)
  diag(sigmas) %*% R %*% diag(sigmas)
}

build_Sigma = function(sigmas,rho,type=NULL){
  if (is.null(type)){
    R = diag(length(sigmas))
    R[R==0] = rho
    outer(sigmas,sigmas) * R
  } else if (type=="ar1") {
    ar1_cov(sigmas, rho)
  } else {
    error("Sigma type not available")
  }
}

sim_Sigma.hat_mv = function(Sigma,nu){
  Sigma.hat = tryCatch({rWishart(n=1,df=nu,Sigma=Sigma)/nu},
                       error=function(e){e})
  if(inherits(Sigma.hat,"error")){
    cat("Error in computing Sigma_nu.hat, probably because p>nu.")
  } else {
    dim(Sigma.hat)=dim(Sigma)
  }
  Sigma.hat
}

sim_theta.hat_mv = function(theta,Sigma){
  p=nrow(Sigma)
  as.vector(rmvnorm(n=1,mean=theta,sigma=Sigma))
}

# compute TOST power
power_TOST_MC_mv = function(alpha, theta, Sigma, nu, delta,
                            B = 10^5, seed = 10^8){
  p=ncol(Sigma)
  set.seed(seed)
  X = rWishart(B, nu, Sigma)/nu
  set.seed(seed*10)
  Z = rmvnorm(n = B, mean = rep(0, p), sigma = Sigma)
  Theta_star=t(as.vector(theta)+t(Z))
  diags=apply(X,3,diag)
  t_val = qt(1 - alpha, df = nu)
  tmp=sqrt(diags)*t_val
  ubs = Theta_star+t(tmp)
  lbs = abs(Theta_star-t(tmp))
  ubs_eval0=(ubs<delta)
  lbs_eval0=(lbs<delta)
  eval=ubs_eval0*lbs_eval0
  list(power_univ=apply(eval,2,mean),power_mult=mean(apply(eval,1,prod)))
}


find_sup_x = function(alpha,Sigma,delta,seed=10^5){
  p=ncol(Sigma)
  if(p<=2){
    method="Brent"
    lower=0;upper=delta
  }else{
    method="Nelder-Mead"
    lower=-Inf;upper=Inf
  }  # We evaluate the size by putting the boundary at every dim
  inds_ = combn(1:p, p-1)
  par0=rep(delta,p-1)
  ## The solution stays sometimes too long at the same place unnecessarily, this is because during optimization
  # the difference between numerical fluctuation of the objective function at the same point
  # are not smaller than tol. --> solution= reltol.
  tmp_=t(apply(inds_,2, function(inds){
    tmp = optim(par=par0,fn=argsup_x,
                inds=inds,alpha=alpha,Sigma=Sigma,delta=delta,seed=seed,
                method=method, lower=lower,upper=upper,
                control=list(reltol=.Machine$double.eps^0.5))
    c(unlist(tmp[1]),unlist(tmp[2]))
    
  }))
  tmp_[,ncol(tmp_)]=round(tmp_[,ncol(tmp_)],4)
  ind_max_size=which.min(tmp_[,ncol(tmp_)])
  inds=inds_[,ind_max_size] #these are the indexes of the dimensions we will let vary, we fix to c the one missing
  thetas=rep(delta,p)
  thetas[inds]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}
# 
argsup_x = function(theta,inds,alpha,Sigma,delta,seed=10^5){
  set.seed(seed)
  p=ncol(Sigma)
  thetas=rep(delta,p)
  thetas[inds]=theta
  # cat(thetas,"\n")
  -power_xTOST_mv(theta = thetas, Sigma = Sigma,
                  delta = delta, alpha=alpha, seed=seed)[1]
}

power_xTOST_mv = function(theta, delta, Sigma, alpha, seed=10^5){
  #delta is what we optimise over
  set.seed(seed)
  mu = rep(0, ncol(Sigma))
  Sig_diag = sqrt(diag(Sigma))
  lb = -delta/Sig_diag-theta/Sig_diag + qnorm(1-alpha)
  ub = delta/Sig_diag-theta/Sig_diag - qnorm(1-alpha)
  if (sum(lb>ub)>0){
    return(0)
  } else {
    tmp_cdf_mvtnorm = mvtnorm::pmvnorm(lower=lb,
                                       upper=ub,
                                       mean=mu,corr=cov2cor(Sigma))[1]
    return(tmp_cdf_mvtnorm)
  }
}


# alphaTOST
obj_fun_alpha_TOST_MC_mv = function(test, alpha, Sigma, nu, delta, theta=NULL, B=10^5, seed=NULL, ...){
  if(is.null(theta)) theta = delta
  if(is.null(seed)){
    size = power_TOST_MC_mv(alpha = test, theta = theta, Sigma = Sigma,
                            nu = nu, delta = delta, B=B)
  }else{
    size = power_TOST_MC_mv(alpha = test, theta = theta, Sigma = Sigma,
                            nu = nu, delta = delta, B=B, seed=seed)
  }
  res = 10000*(size$power_mult - alpha)^2
  res
}

get_alpha_TOST_MC_mv_core = function(alpha, Sigma, nu, delta, theta=NULL, B=10^5, tol = .Machine$double.eps^0.5, seed=NULL, argsup_meth="x", ...){
  if(is.null(seed)) seed=10^8
  if(is.null(theta)){
    if(argsup_meth=="MC"){
      theta = find_sup_MC(alpha,Sigma,nu,delta,B)
    }else if(argsup_meth=="x"){
      theta = find_sup_x(alpha,Sigma,delta)
    }
  }
  obj_func_alpha=obj_fun_alpha_TOST_MC_mv(test = alpha, alpha=alpha, Sigma=Sigma, nu=nu, delta=delta,theta=theta, B=B, seed=seed)
  obj_func_onehalf=obj_fun_alpha_TOST_MC_mv(test = 0.5, alpha=alpha, Sigma=Sigma, nu=nu, delta=delta,theta=theta, B=B, seed=seed)
  if(obj_func_alpha<tol) return(list(minimum=alpha,objective=obj_func_alpha,theta_sup = theta)) #size at a=1 close to alpha, so we consider them equal
  if(obj_func_onehalf<tol) return(list(minimum=alpha,objective=obj_func_onehalf,theta_sup = theta)) #size at a=l close to alpha, so we consider them equal
  out = optimize(obj_fun_alpha_TOST_MC_mv, interval=c(alpha, 0.5), tol=tol,
                 alpha = alpha, Sigma=Sigma, nu = nu,
                 delta = delta, theta = theta, B=B, seed=seed)
  out$theta_sup = theta
  out
}

get_alpha_TOST_MC_mv = function(alpha, Sigma, nu, delta, theta=NULL, B=10^5, tol = .Machine$double.eps^0.5, seed=NULL, max_iter=10, tolpower=NULL, ...){
  if(is.null(tolpower)) tolpower=max(abs(qbinom(c(0.01,0.99),B,alpha)/B-alpha))
  theta_sups = matrix(NA, (max_iter+1), ncol(Sigma))
  if(is.null(theta)) theta_sups[1,]=find_sup_x(alpha,Sigma,delta) else theta_sups[1,]=theta
  alphas = powers = rep(NA, max_iter)
  alphas[1] = alpha
  i=1
  while (i <= max_iter) {
    powers[i] = power_TOST_MC_mv(alpha = alphas[i],
                                 theta = theta_sups[i,], Sigma = Sigma,
                                 nu = nu, delta = delta, B=B)$power_mult
    err_power = abs(powers[i]-alpha)
    if (err_power>tolpower) {
      theta  = theta_sups[i,]
      sol_aTOST = get_alpha_TOST_MC_mv_core(alpha=alpha,
                                            Sigma=Sigma,
                                            nu=nu, delta=delta, theta=theta,
                                            B=B, tol=tol, seed=seed)
      alphas[i+1] = sol_aTOST$min
      theta_sups[i+1,] = find_sup_x(alphas[i+1],Sigma,delta)
      i = i+1
    } else {
      break
    }
  }
  alphas = alphas[!is.na(alphas)]
  powers = powers[!is.na(powers)]
  colnames(theta_sups)=paste0("argsup", 1:ncol(Sigma))
  theta_sups = theta_sups[complete.cases(theta_sups), ]
  if (is.null(nrow(theta_sups))) {
    theta_sup=theta_sups
  } else {
    theta_sup=theta_sups[nrow(theta_sups), ]
  }
  out = list(min=alphas[length(alphas)],
             theta_sup=theta_sup,
             alphas=alphas,
             theta_sups=theta_sups,
             powers=powers,
             iter=i-1,
             err_power=err_power)
  out
}

ci_k = function(theta.hat, Sigma.hat, nu, alpha){
  t_alpha = qt(alpha,df=nu)
  lower <- theta.hat-(-t_alpha)*sqrt(diag(Sigma.hat))
  upper <- theta.hat+(-t_alpha)*sqrt(diag(Sigma.hat))
  # ci = cbind(lower,upper)
  # colnames(ci)=c("lb","ub")
  ci = data.frame(lb=lower,ub=upper)
  ci
}

BE_mv = function(ci,delta){
  return(prod((ci[,1]>-delta)*(ci[,2]<delta)))
}

## compute argsup for theta by MC
get_theta_sup_MC = function(alpha, Sigma, nu, delta, B=10^5, ...){
  # check the sup: which theta has to be on the boundary
  p = ncol(Sigma)
  if (p>1) {
    size_ = rep(0, p)
    for (j in 1:p){
      theta = rep(0,p)
      theta[j] = delta
      size_[j] = power_TOST_MC_mv(alpha = alpha, theta = theta, Sigma = Sigma,
                                  nu = nu, delta = delta, B=B)$power_mult
    }
    if(length(unique(size_))==1) ind=which.max(diag(Sigma)) else ind = which.max(size_)
    theta = rep(0,p)
    theta[ind] = delta
  } else {
    theta=delta
  }
  theta
}


argsup_MC = function(par0,inds,alpha,Sigma,nu,delta,B=10^5,seed=10^8){
  p=ncol(Sigma)
  thetas=rep(delta,p)
  thetas[inds]=par0
  -power_TOST_MC_mv(alpha = alpha, theta = thetas, Sigma = Sigma,
                    nu = nu, delta = delta, B=B, seed=seed)$power_mult
}

find_sup_MC = function(alpha,Sigma,nu,delta,B=10^5,seed=10^8,tolCor=0.01,reltol=.Machine$double.eps^0.5){
  p=ncol(Sigma)
  if (max(cov2cor(Sigma)[upper.tri(Sigma, diag = FALSE)]) < tolCor){
    # correlation quite small: put p-1 thetas at the origin and one on the boundary
    thetas = get_theta_sup_MC(alpha, Sigma, nu, delta, B=B)
    return(thetas)
  }else{
    if(p<=2){
      method="Brent"
      lower=0;upper=delta
    }else{
      method="Nelder-Mead"
      lower=-Inf;upper=Inf
    }
    par0=find_sup_x(alpha,Sigma,delta)
    ind=which.min(abs(par0-delta))
    inds=(1:p)[-ind]
    tmp=optim(par=par0[inds],fn=argsup_MC,
              inds=inds,alpha=alpha,Sigma=Sigma,nu=nu,delta=delta,
              B=B,seed=seed, method=method, lower=lower,upper=upper,
              control=list(reltol=reltol))
    par0[inds]=tmp$par
    return(par0)
  }
}