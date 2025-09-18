# rm(list=ls())
# tikzPlot = 1
main_dir = "simulations/other_plots/"
load_dir = paste0(main_dir, "A6_A7_power_univ_data/")
clusterUse = F

alpha = 0.05
delta = cte = log(1.25)
theta_=seq(0,delta*1.2,by=delta/28)

require(scales)
alp = 1 # transparency
cexsize = 1
sim_all = 1:2
for (sim_setting in sim_all){
  load(paste0(load_dir,"pow_TOST_power_curves.rdata"))
  load(paste0(load_dir,"pow_xTOSTno_power_curves.rdata"))
  load(paste0(load_dir,"pow_xTOSToff_power_curves.rdata"))
  load(paste0(load_dir,"pow_aTOST_power_curves.rdata"))
  nu_all = c(5,10,15,20,25,30,40,50,60)
  if (sim_setting==1){
    nu_ = nu_all[1:3] 
  } else {
    nu_ = nu_all[c(4, 6, 7)]
  }
  # std. error
  sigma_ = c(8, 12, 16)/100
  if (clusterUse==F) {
    if (tikzPlot==1) {
      library(tikzDevice)
      tikz(file=paste0(main_dir, "A6_power_univ_curves_", sim_setting, ".tex"),
           standAlone=TRUE,width=6,height=4.7)
    }
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    color = gg_color_hue(4)
    par(mfrow = c(length(nu_), 3),
        mai = c(0.1, 0.1, 0, 0),
        omi = c(0.6, 0.75, 0.2, 0.05))  # bltr
    for (i in 1:length(sigma_)) {
      for (j in 1:length(nu_)) {
        ind = which(nu_all == nu_[j])
        pow_emp_TOST = 100*pow_TOST_power_curves[ind,i,]
        pow_emp_alphaTOST = 100*pow_aTOST_power_curves[ind,i,]
        pow_emp_xTOST = 100*pow_xTOSTno_power_curves[ind,i,]
        pow_emp_xTOSToff = 100*pow_xTOSToff_power_curves[ind,i,]
        if (j==1) {
          ylim_i = max(100*pow_xTOSTno_power_curves[ind,i,])
        }
        plot(NA, xlim = range(theta_), ylim = c(0, ylim_i), axes = FALSE)
        box()
        grid()
        if (j == 1){
          axis(2, cex.axis = cexsize+0.1)
          mtext(paste0("$ \\sigma_1=", sigma_[i], " $"), cex = 1, line = 2.9, side=2)
        }
        if (i == 1){
          mtext(paste0("$\\nu_2 = ", nu_[j], "$"), cex = 1, line = 0.6, side=3)
        }
        if (i == length(sigma_)){
          axis(1, cex.axis = cexsize+0.1)
        }
        abline(h = 100*alpha, lty = 2, lwd = 1)
        abline(v = delta, lty = 2, lwd = 1)
        lines(theta_, pow_emp_TOST, lwd = 2, col = alpha(color[1], alp), lty = 1)
        lines(theta_, pow_emp_alphaTOST, lwd = 2, col = alpha(color[2], alp), lty = 1)
        lines(theta_, pow_emp_xTOST, lwd = 2, col = alpha(color[3], alp), lty = 1)
        lines(theta_, pow_emp_xTOSToff, lwd = 2, col = alpha(color[4], alp), lty = 1)
        if (i == length(sigma_) && j==2){
          mtext("$\\theta$", cex = cexsize, side=1, line = 3)
        }
        if (j==1 && i==2){
          mtext("Probability of rejecting H$_0$ (\\%)", cex = cexsize, line = 6,
                side=2, adj=0.52)
        }
        if ((i==length(sigma_))&&(j==2)){
          par(xpd=TRUE)
          est_names = c("TOST", "$\\alpha$-TOST", "cTOST", "cTOST*")
          legend("bottomleft",
                 est_names,
                 col = color,
                 lty = 1,
                 lwd = 3, cex = cexsize+0.2,
                 bty = "n",
                 xpd=NA,
                 y.intersp=1,
                 x.intersp=1,
                 xjust=0, yjust=0,
                 text.width=NA,
                 # inset=c(-0.9, -1.1),
                 inset=c(-0.38, -0.6),
                 border = "grey",
                 seg.len=1.75, ncol=4) # horiz = T,
          par(xpd=F)
        }
      }
    }
    if (tikzPlot){
      dev.off()
    } else {
      print(paste0("Figure A", sim_setting+5))
      # invisible(readline(prompt="Press [enter] to continue"))
    }
  }
}
