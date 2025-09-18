library(dplyr)
library(tictoc)
source("aux_fun/ctost_mvt.R")
source("aux_fun/ctost_univ.R")
source("figures/mvt_power_comparison/aux_dTOST.R")
get_params = function(grid, rw){
  grid_rw = grid[rw,]
  sigmas = c(grid_rw$sigma1,grid_rw$sigma2)
  rho = grid_rw$rho
  Sigma = build_Sigma(sigmas=sigmas,rho)
  alpha = grid_rw$alpha
  alpha_t = grid_rw$alpha_t
  return(list(Sigma=Sigma,rho=rho,alpha=alpha,alpha_t=alpha_t))
}
homosk <- c(TRUE, FALSE)
rho <- c(seq(0, 0.95, length.out = 20),0.99)
alpha_t = c(0.01, seq(0.05, 0.45, l=5), 0.5)
alpha <- c(0.05, 0.2)

param_grid_corr <- expand.grid(homosk = homosk, alpha_t = alpha_t, rho = rho, alpha = alpha)

param_grid_corr$sigma1 <- ifelse(param_grid$homosk, 0.075, 0.125)
param_grid_corr$sigma2 <- 0.075

param_grid_corr <- param_grid_corr %>%
  arrange(desc(homosk))

param_grid_corr = param_grid_corr[,-1]

sigma2 <- 0.05
sigma1 <- seq(sigma2, 0.2, by = 0.005)

param_grid_nocorr <- expand.grid(alpha_t = alpha_t, rho = 0, alpha = alpha, sigma1 = sigma1, sigma2 = sigma2)

###
n = nrow(param_grid_corr)+nrow(param_grid_nocorr)

destinationDir <- "power_comparison"
clusterUse <- T

if (clusterUse) {
  if (!dir.exists(destinationDir)) dir.create(destinationDir)
  jobId <- commandArgs(trailingOnly = TRUE)
  jobId <- as.numeric(jobId)
  H <- n  # Total number of simulations
  Rw <- seq(jobId, H, H)  # Distribute jobs across CPUs
} else {
  if(!dir.exists(destinationDir)) dir.create(destinationDir)
  H <- n
  Rw <- 1:H
}

#
nu = 20
cte = log(1.25)
cte_vec = rep(cte, 2)
tolpower=0.0001
B = 1e6

# Start simulations
timing <- rep(NA, length(Rw))
cat("Start -", date(), "\n\n\n")

for(rw in 1:Rw){
  tic()  
  
  res = matrix(NA,nrow=1,ncol=9)
  colnames(res) = c("c1","c2","lamb1","lamb2", "size", "pow_corner", "pow_0", "gamma1", "gamma2")
  rownames(res) = rw
  
  cat("Starting simulation:", which(Rw == rw), "/", length(Rw), "\n")
  if(rw <= nrow(param_grid_corr)) grid = param_grid_corr else grid = param_grid_nocorr
  pars = get_params(grid,rw)
  if(pars$rho==0) lambda = c(cte,0) else lambda = NULL
  #
  dTOST_sol = get_dtost_mvt(alpha=pars$alpha, Sigma=pars$Sigma, delta=cte_vec, theta=lambda,
                            alpha_t = pars$alpha_t, nu=nu, tolpower=tolpower, B=B)
  #
  size = power_TOST_mv(theta=dTOST_sol$theta_sup,
                       delta=dTOST_sol$c_of_0, 
                       Sigma=pars$Sigma,
                       alpha=pars$alpha_t, 
                       nu=nu, 
                       B=B)
  pow_corner = power_TOST_mv(theta=cte_vec,
                             delta=dTOST_sol$c_of_0, 
                             Sigma=pars$Sigma,
                             alpha=pars$alpha_t, 
                             nu=nu, 
                             B=B)
  pow_0 = power_TOST_mv(theta=c(0,0),
                        delta=dTOST_sol$c_of_0, 
                        Sigma=pars$Sigma,
                        alpha=pars$alpha_t, 
                        nu=nu, 
                        B=B)
  gamma1 = power_TOST(alpha=pars$alpha_t, theta=cte,
                      sigma=sqrt(Sigma[1,1]), nu=nu, delta=dTOST_sol$c_of_0[1])
  gamma2 = power_TOST(alpha=pars$alpha_t, theta=cte,
                      sigma=sqrt(Sigma[2,2]), nu=nu, delta=dTOST_sol$c_of_0[2])
  #
  res = c(dTOST_sol$c_of_0, dTOST_sol$theta_sup, size, pow_corner, pow_0, gamma1, gamma2)
  
  # Save results if using cluster
  if (clusterUse) {
    save(res, file = paste0(destinationDir, "/res", rw, ".rdata"))
  }
  
  timing[rw] <- toc()$callback_msg  # Record timing
  cat("Completed simulation:", rw, ";",timing[rw],"\n")
}
































