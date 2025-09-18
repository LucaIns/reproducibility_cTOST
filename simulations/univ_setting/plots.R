# rm(list=ls())
load_dir = "simulations/univ_setting/"
require(tikzDevice)
require(scales)
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
.ep = function(xlim=c(0,1),ylim=c(0,1)) {
  plot(1,1,pch="",axes=FALSE,xlab="",ylab="",main="",
       ylim=ylim,xlim=xlim)
}
load(paste0(load_dir,"pow_TOST_curves_sigma.rdata"))
load(paste0(load_dir,"pow_xTOSTno_curves_sigma.rdata"))
load(paste0(load_dir,"pow_xTOSToff_curves_sigma.rdata"))
load(paste0(load_dir,"pow_aTOST_curves_sigma.rdata"))
# dim is [nu,sigma,theta]
nu_range = as.numeric(gsub("nu=","",dimnames(pow_aTOST_curves_sigma)[[1]]))
sigma_range = as.numeric(gsub("sigma=","",dimnames(pow_aTOST_curves_sigma)[[2]]))
nu_=c(20,40,80)
sigma_ = sigma_range[sigma_range<=0.2]
size_aTOST = pow_aTOST_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),2]
size_TOST = pow_TOST_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),2]
size_cTOST_star = pow_xTOSTno_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),2]
size_cTOST = pow_xTOSToff_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),2]
pow_aTOST = pow_aTOST_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),1]
pow_TOST = pow_TOST_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),1]
pow_cTOST_star = pow_xTOSTno_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),1]
pow_cTOST = pow_xTOSToff_curves_sigma[which(nu_range %in% nu_),which(sigma_range %in% sigma_),1]
if (tikzPlot){
  tikz(paste0(load_dir, "tex/curves_sigma_univ.tex"),
       width = 7, height = 4.5, standAlone = TRUE,
       packages = c(
         "\\usepackage{tikz}",
         "\\usepackage[active,tightpage,psfixbb]{preview}",
         "\\PreviewEnvironment{pgfpicture}",
         "\\setlength\\PreviewBorder{0pt}",
         "\\usepackage{amssymb}",
         "\\usepackage{bm}", "\\usepackage{amsthm}", "\\usepackage{amsbsy}",
         "\\usepackage{amsbsy}",
         "\\usepackage{amsbsy}",
         "\\usepackage{amsfonts}",
         "\\usepackage{amsmath}"
       )
  )
}
my_lwd = 3.5
my_cex_axis = 1.15
B=1e5
sim_error=qbinom(c(0.025,0.975),B,0.05)/B
# Layout
layout(matrix(1:6,ncol=3),width=c(1,1),height=c(1,1))
par(oma=c(2.5,4.7,2.5,0.25))
par(mar=c(1.2,0,0,1))
cols_ = gg_color_hue(4)
sigma1_pow = 0.04
sigma1_size = 0.04
# 1st column 1st row
.ep(c(sigma1_size,max(sigma_)),c(0,0.06))
box()
grid(lty=3,col=alpha("gray",0.4))
legend("bottomleft",legend=c("TOST","$\\alpha$-TOST","cTOST","cTOST*"),
       col=cols_,
       lwd = my_lwd+0.5,
       lty = c(2,3,4,1),
       seg.len = 2.5,
       xpd=NA,
       x.intersp = 0.5,
       bty="n",
       cex=1.2)
u=1
lines(x=c(-16,16),y=c(.05,.05),lty=1,col=alpha(1,0.85))
lines(sigma_, size_TOST[u,], col = cols_[1],lwd=my_lwd,lty=2)
lines(sigma_, size_aTOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, size_cTOST_star[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, size_cTOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
rect(xleft=-16, ybottom=sim_error[1], xright=16, ytop=sim_error[2],border=F,col=alpha(1,0.1))
axis(2,cex.axis=my_cex_axis)
mtext(paste0("$\\nu_2=",nu_[u], "$"),side=3,line=0.8,cex=1.1)
mtext("Empirical size $(\\theta=c_0)$",side=2,line=3.4,cex=1.15,outer=F)
# 1st column 2nd row
.ep(c(sigma1_pow,max(sigma_)),c(0,0.16))
box()
grid(lty=3,col=alpha("gray",0.4))
legend("bottomright",legend=c("$\\alpha$-TOST$-$TOST","cTOST$-$TOST","cTOST*$-$TOST"),
       col=cols_[c(2,3,4)],
       lwd = my_lwd+0.5,
       lty = c(3,4,1),
       seg.len = 2.5,
       xpd=NA,
       x.intersp = 0.5,
       bty="n",
       cex=1.2)
u=1
lines(sigma_, pow_aTOST[u,]-pow_TOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, pow_cTOST_star[u,]-pow_TOST[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, pow_cTOST[u,]-pow_TOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
axis(1,cex.axis=my_cex_axis)
axis(2,cex.axis=my_cex_axis)
mtext("Difference in power $(\\theta=0)$",side=2,line=3.4,cex=1.15,outer=F)
mtext("$\\sigma_1$",side=1,line=2.5,cex=1.2)
# 2nd column 1st row
.ep(c(sigma1_size,max(sigma_)),c(0,0.06))
box()
grid(lty=3,col=alpha("gray",0.4))
u=2
lines(x=c(-16,16),y=c(.05,.05),lty=1,col=alpha(1,0.85))
lines(sigma_, size_TOST[u,], col = cols_[1],lwd=my_lwd,lty=2)
lines(sigma_, size_aTOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, size_cTOST_star[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, size_cTOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
rect(xleft=-16, ybottom=sim_error[1], xright=16, ytop=sim_error[2],border=F,col=alpha(1,0.1))
mtext(paste0("$\\nu_2=",nu_[u], "$"),side=3,line=0.8,cex=1.1)
# 2nd column 2nd row
.ep(c(sigma1_pow,max(sigma_)),c(0,0.16))
box()
grid(lty=3,col=alpha("gray",0.4))
u=2
lines(sigma_, pow_aTOST[u,]-pow_TOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, pow_cTOST_star[u,]-pow_TOST[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, pow_cTOST[u,]-pow_TOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
axis(1,cex.axis=my_cex_axis)
mtext("$\\sigma_1$",side=1,line=2.5,cex=1.2)
# 3rd column 1st row
.ep(c(sigma1_size,max(sigma_)),c(0,0.06))
box()
grid(lty=3,col=alpha("gray",0.4))
u=3
lines(x=c(-16,16),y=c(.05,.05),lty=1,col=alpha(1,0.85))
lines(sigma_, size_TOST[u,], col = cols_[1],lwd=my_lwd,lty=2)
lines(sigma_, size_aTOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, size_cTOST_star[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, size_cTOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
rect(xleft=-16, ybottom=sim_error[1], xright=16, ytop=sim_error[2],border=F,col=alpha(1,0.1))
mtext(paste0("$\\nu_2=",nu_[u], "$"),side=3,line=0.8,cex=1.1)
# 3rd column 2nd row
.ep(c(sigma1_pow,max(sigma_)),c(0,0.16))
box()
grid(lty=3,col=alpha("gray",0.4))
u=3
lines(sigma_, pow_aTOST[u,]-pow_TOST[u,], col = cols_[2],lwd=my_lwd,lty=3)
lines(sigma_, pow_cTOST_star[u,]-pow_TOST[u,], col = cols_[3],lwd=my_lwd,lty=4)
lines(sigma_, pow_cTOST[u,]-pow_TOST[u,], col = cols_[4],lwd=my_lwd,lty=1)
axis(1,cex.axis=my_cex_axis)
mtext("$\\sigma_1$",side=1.5,line=2.5,cex=1.2)
if (tikzPlot){
  dev.off()
} else {
  print("Figure 1")
  # invisible(readline(prompt="Press [enter] to continue"))
}