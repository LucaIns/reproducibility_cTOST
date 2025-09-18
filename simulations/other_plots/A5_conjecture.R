# rm(list=ls())
# tikzPlot = T
main_path = "simulations/other_plots/"
library(dplyr)
library(scales)
B = 1e6
qnt_tol = 0.975
cexsize = 1.5
ylim_diff = c(-0.001, 0.1)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(7)
lwd_lin = 2.5
if (tikzPlot==1) {
  pkgs <- c("\\usepackage{tikz}",
            "\\usepackage[active,tightpage,psfixbb]{preview}",
            "\\PreviewEnvironment{pgfpicture}",
            "\\setlength\\PreviewBorder{0pt}",
            "\\usepackage{amssymb}",
            "\\usepackage{bm}",
            "\\usepackage{amsthm}",
            "\\usepackage{amsbsy}",
            "\\usepackage{amsfonts}")
  library(tikzDevice)
  fig_path = paste0(main_path, "A5_conjecture.tex")
  tikz(file=fig_path,standAlone=TRUE,width=9,height=5,packages=pkgs)
}
layout.matrix <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
layout(mat = layout.matrix, widths = c(1, 1, 1), heights = c(1, 1, 1))
par(mai = c(0.1, 0.3, 0.1, 0), 
    omi = c(0.35, 0.3, 0.21, 1.4))  # bltr
########################
# no-corr
########################
load(file=paste0(main_path, "A5_conjecture_data/param_grid_nocorr.rdata"))
alpha_t = unique(param_grid_nocorr$alpha_t)
alpha = unique(param_grid_nocorr$alpha)
a0 = alpha[1]
dat_heterosk = param_grid_nocorr %>% filter(alpha==a0)
sigma1 = unique(dat_heterosk$sigma1)
sigma2 = unique(dat_heterosk$sigma2)
plot(NA,t="l",xlim=range(sigma1),ylim=range(dat_heterosk$pow_0),
     xlab=NA, ylab=NA, main=NA, axes=F)
grid()
box()
axis(2, cex.axis=cexsize)
axis(1, tick = T, labels = F) 
maintext = paste0("$\\rho=", 0,",\\ \\sigma_{1,2}=",sigma2, "$")
mtext(maintext, cex = cexsize-0.1, line = 1, side=3)
mtext(paste0("Power at $\\bm{\\theta}=\\bm{0}$"), cex = cexsize-0.1, line = 3.75, side=2)
for(a in alpha_t){
  dat = dat_heterosk %>% filter(alpha_t == a)
  lines(sigma1, dat$pow_0,col=cols[which(alpha_t==a)], lwd=lwd_lin)
}
# differences
dat_ref = dat_heterosk %>% filter(alpha_t == 0.5)
dat_ref2 = dat_heterosk %>% filter(alpha_t == 0.45)
diffs = matrix(NA,ncol=length(sigma1),nrow=length(alpha_t)-1)
for(a in alpha_t[-length(alpha_t)]){
  dat = dat_heterosk %>% filter(alpha_t == a) 
  diffs[which(alpha_t==a),]=dat_ref$pow_0-dat$pow_0
} 
tolvar_ref = sapply(dat_ref$pow_0, function(x) 1/B*x*(1-x))
tolvar_2ndref = sapply(dat_ref2$pow_0, function(x) 1/B*x*(1-x))
tol = sapply(sqrt(tolvar_ref+tolvar_2ndref),function(x) c(-1,1)*x*qnorm(qnt_tol))
plot(NA,t="l",xlim=range(sigma1),ylim=ylim_diff,
     xlab=NA, ylab=NA, main=NA, axes=F)
box()
grid()
axis(1, cex.axis=cexsize)
axis(2, tick = T, labels = T, cex.axis=cexsize) 
mtext("$\\sigma_{1,1}$", cex = cexsize-0.1, line = 2.8, side=1)
mtext(paste0("Difference in power at $\\bm{\\theta}=\\bm{0}$"), 
      cex = cexsize-0.1, line = 3.75, side=2)
polygon(x=c(sigma1,rev(sigma1)),
        y=c(tol[2,],rev(tol[1,])),border=F,col=alpha(1,alpha=0.25))
abline(h=0,lty=2)
for(i in 1:(length(alpha_t)-1)){
  lines(sigma1, diffs[i,],col=cols[i], lwd=lwd_lin)
}
########################
# corr
########################
load(file=paste0(main_path, "A5_conjecture_data/param_grid_corr.rdata"))
for (pick_homosk in c(T,F)) {
  head(param_grid_corr)
  alpha_t = unique(param_grid_corr$alpha_t)
  alpha = unique(param_grid_corr$alpha)
  rho = unique(param_grid_corr$rho)
  a0 = alpha[1]
  dat_homosk = param_grid_corr %>% filter(homosk==pick_homosk & alpha==a0)
  sigma1 = unique(dat_homosk$sigma1)
  sigma2 = unique(dat_homosk$sigma2)
  plot(NA,t="l",xlim=range(rho),ylim=range(dat_homosk$pow_0),
       xlab=NA, ylab=NA, main=NA, axes=F)
  grid()
  box()
  axis(2, tick = T, labels = T, cex.axis=cexsize)
  axis(1, tick = T, labels = F) 
  maintext = paste0("$\\sigma_{1,1}=$",sigma1,"$,\ \\sigma_{1,2}=$",sigma2)
  maintext = paste0("$\\sigma_{1,1}=", sigma1,",\\ \\sigma_{1,2}=",sigma2, "$")
  mtext(maintext, cex = cexsize-0.1, line = 1, side=3)
  for(a in alpha_t){
    dat = dat_homosk %>% filter(alpha_t == a) 
    lines(rho, dat$pow_0,col=cols[which(alpha_t==a)], lwd=lwd_lin)
  } 
  if(pick_homosk==F){
    par(xpd=TRUE)
    legend("topright", legend=paste0("$\\alpha_t=",alpha_t, "$"),
           lty=1,col=cols,bty="n",
           xpd=NA, cex = cexsize+0.4,
           x.intersp=0.5,
           xjust=2, lwd=lwd_lin,
           inset=c(-0.65, 0))
    par(xpd=F) 
  }
  # differences
  dat_ref = dat_homosk %>% filter(alpha_t == 0.5)
  dat_ref2 = dat_homosk %>% filter(alpha_t == 0.45)
  diffs = matrix(NA,ncol=length(rho),nrow=length(alpha_t)-1)
  for(a in alpha_t[-length(alpha_t)]){
    dat = dat_homosk %>% filter(alpha_t == a) 
    diffs[which(alpha_t==a),]=dat_ref$pow_0-dat$pow_0
  }
  tolvar_ref = sapply(dat_ref$pow_0, function(x) 1/B*x*(1-x))
  tolvar_2ndref = sapply(dat_ref2$pow_0, function(x) 1/B*x*(1-x))
  tol = sapply(sqrt(tolvar_ref+tolvar_2ndref),function(x) c(-1,1)*x*qnorm(qnt_tol))
  plot(NA,t="l",xlim=range(rho),ylim=ylim_diff,
       xlab=NA, ylab=NA, main=NA, axes=F)
  box()
  grid()
  axis(1, cex.axis=cexsize)
  axis(2, tick = T, labels = F) 
  mtext("$\\rho$", cex = cexsize-0.1, line = 2.8, side=1)
  polygon(x=c(rho,rev(rho)),y=c(tol[2,],rev(tol[1,])),
          border=F,col=alpha(1,alpha=0.25))
  abline(h=0,lty=2)
  for(i in 1:(length(alpha_t)-1)){
    lines(rho, diffs[i,],col=cols[i], lwd=lwd_lin)
  }
}
if (tikzPlot){
  dev.off()
} else {
  print("Figure A5")
  # invisible(readline(prompt="Press [enter] to continue"))
}

