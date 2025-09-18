library(PowerTOST)


sim_theta.hat = function(mu, sigma) rnorm(n=1,mean=mu,sd=sigma)
sim_sigma.hat = function(sigma, nu) sqrt(rchisq(n=1,df=nu)*sigma^2/nu)

power_TOST = function(alpha, theta, sigma, nu, delta, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + delta)/sigma
  mu2 = (theta - delta)/sigma
  R = (delta*sqrt(nu))/(tval*sigma)
  p1 = OwensQ(nu, tval, mu1, 0, R)
  p2 = OwensQ(nu, -tval, mu2, 0, R)
  pw = p2-p1
  pw[pw < 0] = 0
  pw
}

power_xTOST = function(theta, sigma, delta){
  pnorm((theta + delta)/sigma)-pnorm((theta - delta)/sigma)
}

size_xTOST = function(sigma, delta, delta_star, ...){
  power_xTOST(theta = delta, sigma = sigma, delta = delta_star)
}

obj_fun_alpha_TOST = function(test, alpha, sigma, nu, delta, theta=NULL, ...){
  if(is.null(theta)) theta = delta
  size = power_TOST(alpha = test, theta = theta, sigma = sigma,
                    nu = nu, delta = delta)
  (size - alpha)
}

obj_fun_c_of_0 = function(test, alpha, sigma, delta...){
  size = size_xTOST(sigma = sigma, delta = delta, delta_star = test)
  10^6*(size - alpha)
}

get_alpha_TOST = function(alpha, sigma, nu, delta, l=0.5, tol = .Machine$double.eps, ...){
  obj_func_a1=obj_fun_alpha_TOST(test=alpha,alpha=alpha,sigma=sigma, nu=nu, delta=delta)
  obj_func_al=obj_fun_alpha_TOST(test=l,alpha=alpha,sigma=sigma,nu=nu,delta=delta)
  if(obj_func_a1>-tol) return(list(root=alpha,f.root=obj_func_a1)) #size at a=1 close to alpha, so we consider them equal
  if(obj_func_al<tol) return(list(root=l,f.root=obj_func_al)) #size at a=l close to alpha, so we consider them equal
  out = uniroot(obj_fun_alpha_TOST, interval=c(alpha,l),
                alpha = alpha, sigma=sigma,
                nu = nu, delta = delta)
  out
}

get_c_of_0 = function(delta, sigma, alpha, B = 1000, tol = 10^(-8), l=1, optim = "NR"){
  c0=delta
  alpha0=alpha
  # argzero
  if(optim=="uniroot"){
    c_uniroot_ = uniroot(obj_fun_c_of_0, interval=c(10^-8,l),
                         alpha=alpha0, sigma=sigma, delta=c0,tol=.Machine$double.eps)
    size_uniroot = size_xTOST(sigma=sigma,delta=c0,delta_star=c_uniroot_$root)
    out = list(c = c_uniroot_$root, size = size_uniroot)
  }else if(optim=="NR"){
    # Sequence of c's
    m=30
    cte_vect = rep(NA, m)
    # Initial approximation
    c_init = c0 - sigma*qnorm(1 - alpha0)
    if (c_init < 0){
      c_init = c0/2
    }
    cte_vect[1] = c_init
    # Start newton raphson
    for (i in 1:(B-1)){
      delta = (pnorm((c0 + cte_vect[i])/sigma) - pnorm((c0 - cte_vect[i])/sigma) - alpha0)*sigma /
        (dnorm((c0 + cte_vect[i])/sigma) + dnorm((c0 - cte_vect[i])/sigma))
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
      size_NR = size_xTOST(sigma=sigma, delta=c0, delta_star=out$c)
      out = list(c = out$c, size=size_NR, converged = out$converged, iter = out$iter)
    }else{
      size_NR = size_xTOST(sigma=sigma, delta=c0, delta_star=cte_vect[B])
      out = list(c = cte_vect[B], size=size_NR, converged = F, iter = B)
    }
  }
  return(out)
}

ci_BE = function(theta.hat, sigma.hat, nu, alpha){
  t_alpha=qt(p=alpha,df=nu)
  lower <- theta.hat+t_alpha*sigma.hat
  upper <- theta.hat-t_alpha*sigma.hat
  cbind(lower,upper)
}

BE = function(ci, delta){
  return(prod((ci[1]>-delta)*(ci[2]<delta)))
}

BE_TOST = function(theta.hat,sigma.hat,nu,delta,alpha){
  ci = ci_BE(theta.hat=theta.hat,sigma.hat=sigma.hat,alpha=alpha,nu=nu)
  BE(ci,delta)
}

BE_alphaTOST = function(theta.hat,sigma.hat,nu,delta,alpha){
  alpha_star = get_alpha_TOST(alpha=alpha,sigma=sigma,nu=nu,delta=delta)
  ci = ci_BE(theta.hat=theta.hat,sigma.hat=sigma.hat,nu=nu,alpha=alpha_star$root)
  BE_res=BE(ci=ci,delta=delta)
  list(BE_res=BE_res,alpha_star=alpha_star)
}

BE_xTOST = function(theta.hat,sigma.hat,alpha,delta,optim="uniroot", l = 1){
  c_of_0 = get_c_of_0(delta=delta,sigma=sigma.hat,alpha=alpha,optim=optim, l=l)
  BE_res=abs(theta.hat)<c_of_0$c
  list(BE_res=BE_res,c_of_0=c_of_0)
}
