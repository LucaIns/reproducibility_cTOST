rm(list=ls())
tikzPlot = 0
require(scales)
YOUR_PATH = "simulations"
plot_path = paste0(YOUR_PATH, "/tex/mvt_power_sim.tex")
clusterUse = F
alpha = 0.05
delta = cte = log(1.25)
theta_= seq(0, 1.2, length=30)
alp = 1 # transparency
variab_=c("small","medium","large","mix_small","mix_medium","mix_large")
rho_ = c(0, 0.5, 0.9)
if (clusterUse==F) {
  if (tikzPlot==1) {
    library(tikzDevice)
    tikz(file=plot_path,standAlone=TRUE,width=5,height=5.2)
  }
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  color = gg_color_hue(4)[c(1, 2, 3)]
  par(mfrow = c(length(variab_), 3),
      mai = c(0.1, 0.1, 0, 0),
      omi = c(0.6, 0.75, 0.2, 0.05))  # bltr
  for (i in 1:length(variab_)) {
    for (j in 1:3) {
      sim_setting = 1
      load_dir  = paste0(YOUR_PATH, "/mvt_setting/sim_HPC_out/sim_setting_", sim_setting)
      load(file=paste0(load_dir, "/pow_emp_TOST.rdata"))
      load(file=paste0(load_dir, "/pow_emp_xTOST.rdata"))
      if (j==1) {
        jmax = 3
        ylim_i = 100*max(pow_emp_TOST[1,i,jmax,],
                         pow_emp_xTOST[1,i,jmax,]) + 1
      }
      plot(NA, xlim = range(theta_), ylim = c(0, ylim_i), axes = FALSE)
      box()
      grid()
      if (j == 1){
        axis(2)
        mtext(paste0("$\\mathbf{\\Sigma}_1^{(", i, ")} $"), cex = 1, line = 2.2, side=2, las=1)
      }
      if (i == 1){
        mtext(paste0("$ \\rho=", rho_[j], " $"), cex = 1, line = 0.6, side=3)
      }
      if (i == length(variab_)){
        axis(1)
      }
      abline(h = 100*alpha, lty = 2, lwd = 1)
      abline(v = 1, lty = 2, lwd = 1)
      lines(theta_, 100*pow_emp_TOST[1,i,j,], lwd = 2, col = alpha(color[1], alp), lty = 1)
      lines(theta_, 100*pow_emp_xTOST[1,i,j,], lwd = 2, col = alpha(color[3], alp), lty = 1)
      sim_setting = 3
      load_dir  = paste0(YOUR_PATH, "/mvt_setting/sim_HPC_out/sim_setting_", sim_setting)
      load(file=paste0(load_dir, "/pow_emp_TOST.rdata"))
      load(file=paste0(load_dir, "/pow_emp_xTOST.rdata"))
      lines(theta_, 100*pow_emp_TOST[1,i,j,], lwd = 2, col = color[1], lty = 2)
      lines(theta_, 100*pow_emp_xTOST[1,i,j,], lwd = 2, col = color[3], lty = 2)
      cexsize = 1
      if (i == length(variab_) && j==2){
        mtext("$\\kappa$", cex = cexsize, side=1, line = 2.5)
      }
      if (j==1 && i==3){
        mtext("Probability of rejecting H$_0$ (\\%)", cex = cexsize, line = 6,
              side=2, adj=0.9)
      }
      if ((i==length(variab_))&&(j==2)){
        ltys = rep(1, 3)
        par(xpd=TRUE)
        est_names = c("multivariate TOST", "multivariate cTOST")
        myorder = c(1, 3, 2, 4)
        color_red = rep(color[c(1,3)], 2)
        legend("bottomleft",
               c(paste0(est_names, " ($K=2$)"), paste0(est_names, " ($K=4$)"))[myorder],
               col = c(alpha(color_red, alp), color_red)[myorder],
               lty = c(rep(1, 2), rep(2, 2))[myorder],
               lwd = 2.5, cex = cexsize+0.15,
               bty = "n",
               xpd=NA,
               y.intersp=1,
               x.intersp=1,
               xjust=0, yjust=0,
               text.width=NA,
               # inset=c(-0.9, -1.1),
               inset=c(-0.8, -1.2),
               border = "grey",
               seg.len=1.75, ncol=3) # horiz = T,
        par(xpd=F)
      }
    }
  }
  if (tikzPlot){
    dev.off()
  } else {
    print("Figure 3")
    # invisible(readline(prompt="Press [enter] to continue"))
  }
}
