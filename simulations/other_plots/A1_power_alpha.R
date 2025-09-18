# rm(list=ls())
library(PowerTOST)
power_TOST = function(alpha, theta, sigma_nu, nu, delta, delta0, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + delta)/sigma_nu
  mu2 = (theta - delta)/sigma_nu
  R = (delta*sqrt(nu))/(tval*sigma_nu)
  p1 = OwensQ(nu, tval, mu1, 0, R)
  p2 = OwensQ(nu, -tval, mu2, 0, R)
  pw = p2-p1
  pw[pw < 0] = 0
  pw
}
obj_fun_delta_star = function(test, alpha, sigma_nu, nu, delta, ...){
  size = power_TOST(alpha = alpha, theta = delta, sigma_nu = sigma_nu,
                    nu = nu, delta = test)
  (size - 0.05) # hard-coded alpha0
}
get_delta_star = function(delta, alpha, sigma_nu, nu, l=100, alpha0, tol = 1e-16, ...){
  out = uniroot(obj_fun_delta_star,interval=c(1e-6,l),
                alpha = alpha, sigma_nu = sigma_nu,
                nu = nu, delta = delta, alpha0=alpha0, tol = tol)
  out
}
power_normal = function(alpha,theta,c,sigma_nu){
  zalpha = qnorm(1-alpha)
  b = (theta/sigma_nu +c(-1,1)*(c/sigma_nu - zalpha))
  pnorm(b[2])-pnorm(b[1])
}
obj_fun_delta_star_normal = function(test,alpha,c0,sigma_nu){
  size = power_normal(alpha=alpha,theta=c0,c=test,sigma_nu=sigma_nu)
  (size - 0.05)
}
get_delta_star_normal = function(alpha, sigma_nu, c0=c0, l=100, tol = 1e-16, ...){
  out = uniroot(obj_fun_delta_star_normal,interval=c(1e-6,l),
                alpha = alpha, sigma_nu = sigma_nu, c0 = c0, tol = tol)
  out
}
nu=15
sigma_nu=0.1
c0 = log(1.25)
alpha0=0.05
alpha_ = seq(0.001,1/2,l=100)
cstar=sapply(alpha_,function(x) get_delta_star(c0,x,sigma_nu,nu)$root)
cstarn = sapply(alpha_,function(x) get_delta_star_normal(x,sigma_nu,c0)$root)
power = sapply(1:length(alpha_),function(x) power_TOST(alpha_[x],0,sigma_nu,nu,cstar[x],c0))
powern = sapply(1:length(alpha_),function(x) power_normal(alpha_[x],0,cstarn[x],sigma_nu))
# install.packages("devtools")
# devtools::install_github("stephaneguerrier/cTOST")
require(cTOST)
atost = ctost(theta = 0, sigma = sigma_nu^2, nu = nu, 
              delta = c0, method = "alpha")
atost$corrected_alpha
library(RColorBrewer)
color = brewer.pal(8, "Dark2")[3:8]
library(tikzDevice)
if (tikzPlot){
  mypath_output = "./simulations/other_plots/A1_power_alpha.tex"
  tikz(file=mypath_output,standAlone=TRUE,width=6,height=4)
}
par(mfrow = c(1, 1), 
    mai = c(1.2, 0.75, 0.1, 0.1))  # bltr
plot(NA,xlim=range(alpha_),ylim=range(c(powern,power)),main="",
     xaxt="n", yaxt="n", ylab="", xlab="") # nu=15, sigma=0.1
axis(2,cex.axis=1.2)
mtext("$\\alpha$", side=1, line=2.7, cex=1.5)
axis(1,cex.axis=1.2)
mtext("Power at $\\theta=0$", side=2, line=2.7, cex=1.5)
lines(alpha_,powern,t="l",col=color[1], lwd=4)
lines(alpha_,power,t="l",col=color[2], lwd=4)
abline(v=0.05, col=color[3], lty=2, lwd=4)
abline(v=atost$corrected_alpha, col=color[4], lty=2, lwd=4)
abline(v=0.5, col=color[5], lty=2, lwd=4)
grid()
legend("bottom", legend=c("Known $\\sigma_1$", "Unknown $\\sigma_1$"),
       col=c(color[1], color[2]), lty=c(1,1), lwd=4, cex=1.25,
       horiz=TRUE, inset=c(0,-0.37), xpd=TRUE,  bty = "n",  x.intersp=0.5)
legend("bottom", legend=c("$\\alpha=0.05$", 
                          paste0("$\\alpha \\approx $", " ", round(atost$corrected_alpha, 3)), 
                          "$\\alpha=0.5$"),
       col=c(color[3], color[4], color[5]), lty=c(2,2,2), lwd=4, cex=1.25,
       horiz=TRUE, inset=c(0,-0.47), xpd=TRUE,  bty = "n", x.intersp=0.5, 
       xjust=0, yjust=0,
       text.width=c(0.145,0.155,0.12)*0.50)
if (tikzPlot){
  dev.off()
} else {
  print("Figure A1")
  # invisible(readline(prompt="Press [enter] to continue"))
}