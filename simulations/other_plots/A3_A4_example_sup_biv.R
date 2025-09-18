# rm(list=ls())
plot_path = "simulations/other_plots"
library(tikzDevice)
library(scales)
size1 = function(c1, theta1, sigma1){
  (pnorm((theta1 + c1)/sigma1) - pnorm((theta1 - c1)/sigma1))
}
obj1 = function(c1, theta1, sigma1, alpha0){
  (size1(c1, theta1, sigma1) - alpha0)^2
}
size2 = function(c1, c2, theta1, theta2, sigma1, sigma2){
  (pnorm((theta1 + c1)/sigma1) - pnorm((theta1 - c1)/sigma1))*(pnorm((theta2 + c2)/sigma2) - pnorm((theta2 - c2)/sigma2))
}
obj2 = function(c1, c2, c0, sigma1, sigma2, alpha0){
  s1 = size2(c1 = c1, c2 = c2, theta1 = c0, theta2 = 0, sigma1 = sigma1, sigma2 = sigma2) 
  s2 = size2(c1 = c1, c2 = c2, theta1 = 0, theta2 = c0, sigma1 = sigma1, sigma2 = sigma2) 
  (max(s1, s2) - alpha0)^2
}
obj_sup = function(x, c1, c2, sigma1, sigma2){
  if (max(abs(x[1]), abs(x[2])) < (log(1.25) - 0.01)){
    10
  }else{
    -size2(c1 = c1, c2 = c2, theta1 = x[1], theta2 = x[2], sigma1 = sigma1, sigma2 = sigma2)
  }
}
check_sup = function(c0, c1, c2, sigma1, sigma2){
  c01 = seq(from = c0, to = -c0, length.out = 100)
  c02 = rep(c0, 100)
  tests1 = cbind(c01, c02)
  tests2 = cbind(c01, -c02)
  tests3 = cbind(c02, c01)
  tests4 = cbind(-c02, c01)
  tests = rbind(tests1, tests2, tests3, tests4)
  m = nrow(tests)
  res = rep(NA, m)
  for (i in 1:m){
    res[i] = size2(c1 = c1, c2 = c2, theta1 = tests[i,1], theta2 = tests[i,2], sigma1 = sigma1, sigma2 = sigma2)
  }
  max(res)
}
c0 = log(1.25)
alpha0 = 0.05
sigma1_all = c(0.1, 0.1, 0.1)
sigma2_all = c(0.1, 0.12, 0.14)
B = 1000
cexsize = 1.5
linedist = 3.2
point_cex = 2.2
lwd_lin = 2.5
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
######################################
# Figure A3
######################################
cols_lin = gg_color_hue(4)
if (tikzPlot){
  tikz(paste0(plot_path, "/A3_example_sup_biv_curves.tex"), width = 8, height = 8, standAlone = TRUE,
       packages = c("\\usepackage{tikz}",
                    "\\usepackage[active,tightpage,psfixbb]{preview}",
                    "\\PreviewEnvironment{pgfpicture}",
                    "\\setlength\\PreviewBorder{0pt}",
                    "\\usepackage{amssymb}",
                    "\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}"
                    ,"\\usepackage{amsbsy}"
                    ,"\\usepackage{amsbsy}"
                    ,"\\usepackage{amsfonts}"))
}
par(mfrow=c(3,3),
    mar=c(1.5,1,1,1)*0.5,
    oma=c(3.75,4.4,3.5,0))
layout.matrix = matrix(1:9, nrow = 3, ncol = 3)
layout(mat = layout.matrix)
for (h in 1:3){
  sigma1 = sigma1_all[h]
  sigma2 = sigma2_all[h]
  c11 = optimize(obj1, c(0, 1), theta1 = c0, sigma1 = sigma1, alpha0 = alpha0)$minimum
  c21 = optimize(obj1, c(0, 1), theta1 = c0, sigma1 = sigma2, alpha0 = alpha0)$minimum
  size1(c11, theta1 = c0, sigma1 = sigma1)
  size1(c21, theta1 = c0, sigma1 = sigma2)
  c22 = seq(from = 0.025, to = 0.145, length.out = B)
  c12 = which_sup = size = power = marginal_size1 = marginal_size2 = rep(NA, B)
  for (i in 1:B){
    c12[i] = optimize(obj2, c(0, 1), c2 = c22[i], c0 = c0, sigma1 = sigma1, sigma2 = sigma2, alpha0 = alpha0)$minimum
    inter = check_sup(c0 = c0, c1 = c12[i], c2 = c22[i], sigma1 = sigma1, sigma2 = sigma2)
    if (abs(inter[1] - 0.05) > 0.01){
      print(i)
    }
    # Size
    size_inter = c(size2(c1 = c12[i], c2 = c22[i], theta1 = c0, theta2 = 0, sigma1 = sigma1, sigma2 = sigma2),
                   size2(c1 = c12[i], c2 = c22[i], theta1 = 0, theta2 = c0, sigma1 = sigma1, sigma2 = sigma2))
    which_sup[i] = which.max(size_inter)
    size[i] = size_inter[which_sup[i]]
    power[i] = size2(c1 = c12[i], c2 = c22[i], theta1 = 0, theta2 = 0, sigma1 = sigma1, sigma2 = sigma2)
    marginal_size1[i] = size1(c1 = c12[i], theta1 = c0, sigma1 = sigma1)
    marginal_size2[i] = size1(c1 = c22[i], theta1 = c0, sigma1 = sigma2)
  }
  plot(c22, c12, type = "l", lwd=lwd_lin,
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = 1)
  if (h==1) {
    axis(2, cex.axis=cexsize)
    mtext("$c_1$", side = 2, line = linedist, cex=cexsize)
  }
  rect(-1, -1, c22[1 + which.max(diff(which_sup))], 2, col = "grey98", border = NA)
  rect(c22[1 + which.max(diff(which_sup))], -1, 1, 2, col = "grey92", border = NA)
  lines(c22, c12, lwd=lwd_lin)
  grid()
  box()
  points(c22[which.min(abs(marginal_size1 - marginal_size2))], c12[which.min(abs(marginal_size1 - marginal_size2))], 
         pch = 16, cex = point_cex, col=alpha("orange", 1))
  points(c22[which.max(power)], c12[which.max(power)], 
         pch = 1, cex = point_cex, lwd=2, col="blue")  
  tit = paste0("$\\sigma_{1,1}=", sigma1, "$", ", $\\sigma_{1,2}=", sigma2, "$")
  mtext(tit, side = 3, line = 1, cex=cexsize)
  ###############
  # FIG 2
  ###############
  plot(c22, power, type = "l", col = 1, ylim = c(0.1,0.4),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", lwd=lwd_lin)
  if (h==1) {
    axis(2, cex.axis=cexsize)
    mtext("Power at (0,0)", side = 2, line = linedist, cex=cexsize)
  }
  rect(-1, -1, c22[1 + which.max(diff(which_sup))], 2, col = "grey98", border = NA)
  rect(c22[1 + which.max(diff(which_sup))], -1, 1, 2, col = "grey92", border = NA)
  lines(c22, power, lwd=lwd_lin)
  grid()
  box()
  grid()
  points(c22[which.min(abs(marginal_size1 - marginal_size2))], 
         power[which.min(abs(marginal_size1 - marginal_size2))], 
         pch = 16, cex = point_cex, col=alpha("orange", 1))
  points(c22[which.max(power)], power[which.max(power)], 
         pch = 1, cex = point_cex, lwd=2, col="blue")
  plot(c22, power, type = "l", col = 1, ylim = c(0,0.3),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", lwd=lwd_lin)
  mtext("$c_2$", side = 1, line = linedist, cex=cexsize)
  if (h==1) {
    axis(2, cex.axis=cexsize)
    mtext("Marginal size", side = 2, line = linedist, cex=cexsize)
  }
  axis(1, cex.axis=cexsize)
  rect(-1, -1, c22[1 + which.max(diff(which_sup))], 2, col = "grey98", border = NA)
  rect(c22[1 + which.max(diff(which_sup))], -1, 1, 2, col = "grey92", border = NA)
  lines(c22, marginal_size1, col = cols_lin[1], lwd=lwd_lin)
  lines(c22, marginal_size2, col = cols_lin[2], lwd=lwd_lin)
  abline(h = 0.05, lty = 2)
  box()
  grid()
  points(c22[which.min(abs(marginal_size1 - marginal_size2))], 
         marginal_size1[which.min(abs(marginal_size1 - marginal_size2))], 
         pch = 16, cex = point_cex, col=alpha("orange", 1))
  points(c22[which.max(power)], marginal_size1[which.max(power)], 
         pch = 1, cex = point_cex, lwd=2, col="blue")
  points(c22[which.max(power)], marginal_size2[which.max(power)], 
         pch = 1, cex = point_cex, lwd=2, col="blue")
}
if (tikzPlot){
  dev.off()
} else {
  print("Figure A3")
  # invisible(readline(prompt="Press [enter] to continue"))
}
######################################
# Figure A4
######################################
if (tikzPlot){
  tikz(paste0(plot_path, "/A4_example_sup_biv_rej_reg.tex"), width = 8, height = 4, standAlone = TRUE,
       packages = c("\\usepackage{tikz}",
                    "\\usepackage[active,tightpage,psfixbb]{preview}",
                    "\\PreviewEnvironment{pgfpicture}",
                    "\\setlength\\PreviewBorder{0pt}",
                    "\\usepackage{amssymb}",
                    "\\usepackage{bm}","\\usepackage{amsthm}","\\usepackage{amsbsy}"
                    ,"\\usepackage{amsbsy}"
                    ,"\\usepackage{amsbsy}"
                    ,"\\usepackage{amsfonts}"))
}
par(mfrow=c(1,2),
    mar=c(1,1,1,1),
    oma=c(4,0,1.2,0))
layout(mat = matrix(1:2, ncol=2))
for (h in 2:3){
  sigma1 = sigma1_all[h]
  sigma2 = sigma2_all[h]
  c11 = optimize(obj1, c(0, 1), theta1 = c0, sigma1 = sigma1, alpha0 = alpha0)$minimum
  c21 = optimize(obj1, c(0, 1), theta1 = c0, sigma1 = sigma2, alpha0 = alpha0)$minimum
  size1(c11, theta1 = c0, sigma1 = sigma1)
  size1(c21, theta1 = c0, sigma1 = sigma2)
  c22 = seq(from = 0.01, to = 0.2, length.out = B)
  c12 = which_sup = size = power = marginal_size1 = marginal_size2 = rep(NA, B)
  for (i in 1:B){
    c12[i] = optimize(obj2, c(0, 1), c2 = c22[i], c0 = c0, sigma1 = sigma1, sigma2 = sigma2, alpha0 = alpha0)$minimum
    inter = check_sup(c0 = c0, c1 = c12[i], c2 = c22[i], sigma1 = sigma1, sigma2 = sigma2)
    if (abs(inter[1] - 0.05) > 0.01){
      print(i)
    }
    size_inter = c(size2(c1 = c12[i], c2 = c22[i], theta1 = c0, theta2 = 0, sigma1 = sigma1, sigma2 = sigma2),
                   size2(c1 = c12[i], c2 = c22[i], theta1 = 0, theta2 = c0, sigma1 = sigma1, sigma2 = sigma2))
    which_sup[i] = which.max(size_inter)
    size[i] = size_inter[which_sup[i]]
    power[i] = size2(c1 = c12[i], c2 = c22[i], theta1 = 0, theta2 = 0, sigma1 = sigma1, sigma2 = sigma2)
    marginal_size1[i] = size1(c1 = c12[i], theta1 = c0, sigma1 = sigma1)
    marginal_size2[i] = size1(c1 = c22[i], theta1 = c0, sigma1 = sigma2)
  }
  ###############
  # FIG 
  ############### 
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
  cols = gg_color_hue(3)
  cols = cols[-2]
  reds = t_col(cols[1], c(0, 20, 40, 85))
  blues = t_col(cols[2], c(0, 20, 40, 80))
  col1 = "grey"
  grey = t_col(col1, 80)
  grey2 = t_col(col1, 30)
  sig = 0.05
  delta = sig*qnorm(1-0.05)
  c_val = log(1.25)
  alphas = c(0.2, 0.05, 0.01)
  K = length(alphas)
  eps = 0.0075
  plot(NA, xlim = c(-0.25, 0.25), ylim = c(-0.25, 0.25),
       axes = FALSE, xlab = "", ylab = "")
  rect(-c_val, -c_val, c_val, c_val, col = grey, border = grey2)
  require(scales)
  delta1 = sigma1*qnorm(1-0.05)
  delta2 = sigma2*qnorm(1-0.05)
  rect(-c_val + delta1, -c_val + delta2, c_val - delta1, c_val - delta2, 
       col = alpha("purple", 0.25), border = "purple")
  # our solution
  c2star = c22[which.min(abs(marginal_size1 - marginal_size2))]
  c1star = c12[which.min(abs(marginal_size1 - marginal_size2))]
  rect(-c1star, -c2star, c1star, c2star, col = alpha("orange", 0.5), border = "orange", lty=1, density=30, angle=-45)
  cols = c("purple", "orange", "blue")
  # max power
  c2star = c22[which.max(power)]
  c1star = c12[which.max(power)]
  rect(-c1star, -c2star, c1star, c2star, col = alpha("blue", 0.5), border = "blue", lty=1, density=30)
  cexsize = 1.2
  lines(c(-1, 1), c(0, 0))
  lines(c(0, 0), c(-1, 1))
  lines(c(c_val, c_val), c(-eps, eps))
  lines(c(-c_val, -c_val), c(-eps, eps))
  lines(c(-eps, eps), c(c_val, c_val))
  lines(c(-eps, eps), c(-c_val, -c_val))
  text(c_val, -2*eps, "$c_0$", cex=cexsize)
  text(-c_val, -2*eps, "$-c_0$", cex=cexsize)
  text(-2.5*eps, c_val+eps/2, "$c_0$", cex=cexsize)
  text(-2.85*eps, -c_val-eps/2,  "$-c_0$", cex=cexsize)
  text(1.9*eps, -1.9*eps, "$0$", cex=cexsize)
  text(0.28, -2.5*eps, "$\\theta_1 $", cex = cexsize+0.2, xpd=NA)
  text(-2.5*eps, -0.29, "$\\theta_2 $", cex = cexsize+0.2, xpd=NA)
  tit = paste0("$\\sigma_{1,1}=", sigma1, "$, ", 
               "$\\sigma_{1,2}=", sigma2, "$")
  mtext(tit, side = 3, line = 1, cex=cexsize+0.2)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), 
    mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(x="bottom", inset = 0,
       legend=c("TOST", "cTOST", "other"),
       fill=c(alpha(cols, 0.6)),
       cex=cexsize+0.2,
       horiz=TRUE,
       bty = "n", border = NA,
       title = "Rejection regions")
if (tikzPlot){
  dev.off()
} else {
  print("Figure A4")
  # invisible(readline(prompt="Press [enter] to continue"))
}