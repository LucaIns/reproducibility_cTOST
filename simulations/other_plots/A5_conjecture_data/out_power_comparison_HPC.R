rm(list=ls())
# library(dplyr)
homosk <- c(TRUE, FALSE)
rho <- c(seq(0, 0.95, length.out = 20),0.99)
alpha_t = c(0.01, seq(0.05, 0.45, l=5), 0.5)
alpha <- c(0.05, 0.2)
#
param_grid_corr <- expand.grid(homosk = homosk, alpha_t = alpha_t, rho = rho, alpha = alpha)
#
param_grid_corr$sigma1 <- ifelse(param_grid_corr$homosk, 0.075, 0.125)
param_grid_corr$sigma2 <- 0.075
#
param_grid_corr <- param_grid_corr %>%
  arrange(desc(homosk))
#
sigma2 <- 0.05
sigma1 <- seq(sigma2, 0.2, by = 0.005)
#
param_grid_nocorr <- expand.grid(alpha_t = alpha_t, rho = 0, alpha = alpha, sigma1 = sigma1, sigma2 = sigma2)
#
destinationDir <- "power_comparison"
incomplete_out = c()
#
n = nrow(param_grid_corr)+nrow(param_grid_nocorr)
#
param_grid_corr[c("c1","c2","lamb1","lamb2", "size", "pow_corner", "pow_0", "gamma1", "gamma2")] <- NA
param_grid_nocorr[c("c1","c2","lamb1","lamb2", "size", "pow_corner", "pow_0", "gamma1", "gamma2")] <- NA
#
for(rw in 1:n){
  print(rw)
  res = tryCatch({get(load(paste0(destinationDir, "/res", rw, ".rdata")))},error=function(e){e})
  if(inherits(res,"error")){
    incomplete_out = c(incomplete_out,rw)
    next
  }
  if(rw <= nrow(param_grid_corr)){
    param_grid_corr[rw,7:ncol(param_grid_corr)] <- res
  } else {
    rrw=rw-nrow(param_grid_corr)
    param_grid_nocorr[rrw,6:ncol(param_grid_nocorr)] <- res
  } 
}

save(param_grid_corr,file="param_grid_corr.rdata")
save(param_grid_nocorr,file="param_grid_nocorr.rdata")



