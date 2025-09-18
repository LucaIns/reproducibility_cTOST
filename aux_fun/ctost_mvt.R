require(mvtnorm)

# Functions to generate Sigma and Sigma hat
ar1_cor <- function(p, rho) {
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

# compute TOST power
power_TOST_MC_mv = function(alpha, theta, Sigma, nu, delta,
                            B = 10^5, seed = 13847){
  p=ncol(Sigma)
  set.seed(seed)
  X = rWishart(B, nu, Sigma)/nu
  Z = rmvnorm(n = B, mean = rep(0, p), sigma = Sigma)
  Theta_star=t(as.vector(theta)+t(Z))
  diags=apply(X,3,diag)
  t_val = qt(1 - alpha, df = nu)
  tmp=sqrt(diags)*t_val
  ubs = Theta_star+t(tmp)
  lbs = abs(Theta_star-t(tmp))
  ubs_eval0=(t(ubs)<delta)
  lbs_eval0=(t(lbs)<delta)
  eval=t(ubs_eval0*lbs_eval0)
  list(power_univ=apply(eval,2,mean),power_mult=mean(apply(eval,1,prod)))
}

# sup tost
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
  par0=rep(log(1.25),p-1)
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
  thetas=rep(log(1.25),p)
  thetas[inds]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}

argsup_x = function(theta,inds,alpha,Sigma,delta,seed=10^5){
  set.seed(seed)
  p=ncol(Sigma)
  thetas=rep(log(1.25),p)
  thetas[inds]=theta
  -power_cTOST_mv(theta = thetas, Sigma = Sigma,
                  delta = delta, alpha=alpha, seed=seed)[1]
}

# power ctost
power_cTOST_mv = function(theta, delta, Sigma, alpha=1/2, seed=10^5){
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

# sup ctost
find_sup_ctost = function(Sigma,delta,c_of_0){
  p=ncol(Sigma)
  if(p<=2){
    method="Brent"
  }else{
    method="Nelder-Mead"
    # these are the only possible values for lower and upper
    lower=-Inf;upper=Inf
  }
  # dimss2free yields all the combinations for the dimensions that we would let free.
  dims2free_ = combn(1:p, p-1)

  ## The solution stays sometimes too long at the same place unnecessarily, this is because during optimization
  # the difference between numerical fluctuation of the objective function at the same point
  # are not smaller than tol. --> solution= reltol.
  #
  # we optimize for each combination of free dimensions. tmp_ outputs the parameter values
  # of the argsup's free dimension, and the value of the objective function (negative size)
  tmp_=t(apply(dims2free_,2, function(dims2free){
    if(p<=2){
      lower=-delta[dims2free];upper=delta[dims2free]
    }
    tmp = optim(par=rep(0,p-1),fn=objfun4sup_ctost,
                dims2free=dims2free,Sigma=Sigma,delta=delta, c_of_0=c_of_0,
                method=method, lower=lower,upper=upper,
                control=list(reltol=.Machine$double.eps^0.5))
    c(unlist(tmp[1]),unlist(tmp[2])) # the argsups fluctuate a bit but the objective function value is the very close each time
  }))
  # we look from among all the combinations, the optimal yielding the min negative size
  # tmp_[,ncol(tmp_)]=round(tmp_[,ncol(tmp_)],4)
  ind_max_size=which.min(tmp_[,ncol(tmp_)])
  dims2free=dims2free_[,ind_max_size]
  thetas=delta
  thetas[dims2free]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}

objfun4sup_ctost = function(theta,dims2free,Sigma,delta,c_of_0){
  p=ncol(Sigma)
  # we define the potential argsup as a vector where we fix one dimension at
  # its corresponding delta, and let the other dimensions be fed by the values
  # considered by the optimisation
  thetas=delta # this is only a useful trick to get directly the dimension to fix at its corresponding delta
  thetas[dims2free]=theta # then we replace the dimensions we let free (to optimize) by the values considered by the optimization
  # cat(thetas,"\n")
  10^3*(-power_cTOST_mv(theta = thetas, Sigma = Sigma,
                        delta = c_of_0)[1])
}

# multiv cTOST
get_starting_values = function(Sigma, cte, alpha){
  p = ncol(Sigma)
  c_start = rep(NA, p)
  for (i in 1:p){
    c_start[i] = get_c_of_0(cte[i], sqrt(Sigma[i,i]), alpha)$c
  }
  if (any(c_start>cte)) warning("c_start>cte")
  c_start
}

obj_fun_cTOST_constr_size_mv = function(test, Sigma, cte, delta, alpha){
  
  p = ncol(Sigma)
  delta_vec = rep(max(cte), p) # symmetric margins (extract largest lambda)
  eval_f = get_starting_values(Sigma, delta_vec, test)
  tmp_cdf = power_cTOST_mv(cte, eval_f, Sigma)
  10^5*(tmp_cdf - alpha)^2
}


get_c_of_0_mvt = function(start, Sigma, cte, alpha, B1=100, B2=2000,
                          interval_size1=c(0,1), interval_size2=c(0,5)){
  p = ncol(Sigma)
  out = optim(1, obj_fun_cTOST_constr_size_mv,
              alpha = alpha, Sigma = Sigma, delta=start,
              cte = cte, method = "Brent", lower = alpha,
              upper =  1/2)$par
  out2 = list()
  out2$par_mult = out
  delta_vec = rep(max(cte), ncol(Sigma))
  out2$par = get_starting_values(Sigma, delta_vec, out)
  out2
}


get_ctost_mvt = function(alpha, Sigma, delta,
                         theta=NULL, tol = .Machine$double.eps^0.5,
                         seed=NULL, max_iter=10, tolpower=NULL, ...){

  B = 10^3
  if(is.null(tolpower)) tolpower=max(abs(qbinom(c(0.01,0.99),B,alpha)/B-alpha))
  theta_sups = matrix(NA, (max_iter+1), ncol(Sigma))
  if(is.null(theta)) theta_sups[1,]=find_sup_ctost(Sigma,delta,delta) else theta_sups[1,]=theta
  powers = rep(NA, max_iter)
  c_of_0 = matrix(NA, max_iter+1, ncol(Sigma))
  c_of_0[1,] = delta

  i=1
  while (i <= max_iter) {
    powers[i] = power_cTOST_mv(theta=theta_sups[i,], delta=c_of_0[i,], Sigma=Sigma)
    err_power = abs(powers[i]-alpha)
    if (err_power>tolpower) {
      print(paste0("iteration mvt xTOST: ", i))
      theta  = theta_sups[i,]
      sol_xTOST = get_c_of_0_mvt(start=get_starting_values(Sigma=Sigma, cte=delta, alpha=alpha),
                                 Sigma=Sigma,
                                 cte=theta,
                                 alpha=alpha)
      c_of_0[i+1, ] = sol_xTOST$par
      theta_sups[i+1,] = find_sup_ctost(Sigma, delta, c_of_0[i+1, ])
      i = i+1
    } else {
      break
    }
  }
  c_of_0 = matrix(c_of_0[!is.na(c_of_0[,1]),], ncol=ncol(Sigma))
  powers = powers[!is.na(powers)]
  colnames(theta_sups)=paste0("argsup", 1:ncol(Sigma))
  theta_sups = theta_sups[complete.cases(theta_sups), ]
  if (is.null(nrow(theta_sups))) {
    theta_sup=theta_sups
  } else {
    theta_sup=theta_sups[nrow(theta_sups), ]
  }
  out = list(c_of_0=c_of_0[nrow(c_of_0),],
             theta_sup=theta_sup,
             c_of_0_seq=c_of_0,
             theta_sups=theta_sups,
             powers=powers,
             iter=i-1,
             err_power=err_power)
  out
}

