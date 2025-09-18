# rm(list=ls())
par(mfrow=c(1,1))
main_dir = "simulations/other_plots/"
load_dir = paste0(main_dir, "A8_A11_size_heatmap_data/")
library(scales)
library(tikzDevice)
load(paste0(load_dir,"pow_TOST_tiles.rdata"))
load(paste0(load_dir,"pow_xTOSTno_tiles.rdata"))
load(paste0(load_dir,"pow_xTOSToff_tiles.rdata"))
load(paste0(load_dir,"pow_aTOST_tiles.rdata"))
# select nu from 20 to 80
nu_range=16:76
# compute the quantiles and have common quantiles for each method
dats = rbind(pow_TOST_tiles[nu_range,,2],
             pow_xTOSTno_tiles[nu_range,,2],
             pow_xTOSToff_tiles[nu_range,,2],
             pow_aTOST_tiles[nu_range,,2])
# dat is what the tiles will be computed for
# size is 2, power is 1 in last dimension
for (d_typ in 1:4){
  if (d_typ == 1){
    dat = pow_TOST_tiles[nu_range,,2]
    label="fig_size_tost_tiles.tex"
  } else if (d_typ == 2){
    dat = pow_aTOST_tiles[nu_range,,2]
    label="fig_size_atost_tiles.tex"
  } else if (d_typ == 3){
    dat = pow_xTOSTno_tiles[nu_range,,2]
    label="fig_size_ctost_tiles.tex"
  } else if (d_typ == 4){
    dat = pow_xTOSToff_tiles[nu_range,,2]
    label="fig_size_ctost_ref_tiles.tex"
  }
  alpha = 0.05
  delta = log(1.25)
  gamma_=c(1,1.25)
  theta_=log(gamma_)
  nu_ = seq(20,80)
  sigma_ = seq(0.01,0.3,by=0.005)
  nu = nu_
  sig = sigma_
  n1 = length(nu)
  n2 = length(sig)
  alpha = 0.05
  lower = qbinom(.005,100000,alpha)/100000
  upper = qbinom(.995,100000,alpha)/100000
  q1 = quantile(dats[dats <= lower], na.rm = T, probs = 1:6/6)
  q2 = quantile(dats[dats >= lower & dats <= upper], na.rm = T, probs = 1:3/3)
  q3 = quantile(dats[dats > upper], na.rm = T, probs = 1:2/2)
  c(q1,q2,q3)
  diff(sapply(as.vector(c(q1,q2,q3)),function(x) mean(as.vector(dat)<=x)))
  r1 = colorRampPalette(c("#154360", "#85c1e9"))
  col1 = r1(6)
  r2 = colorRampPalette(c("#d6eaf8", "#fadbd8"))
  col2 = r2(3)
  r3 = colorRampPalette(c("#f1948a", "darkred"))
  col3 = r3(2)
  index = matrix(NA, n1, n2)
  for (i in 1:n1){
    for (j in 1:n2){
      val = dat[i,j]
      if (is.na(val)){
      }else{
        if (val < lower){
          index[i,j] = col1[min(which(val<=q1))]
        }else{
          if(val <= upper){
            index[i,j] = col2[min(which(val<=q2))]
          }else{
            index[i,j] = col3[min(which(val<=q3))]
          }
        }
      }
    }
  }
  if (tikzPlot){
    library(tikzDevice)
    tikz(paste0(main_dir, "A8_A11_size_heatmap_", label), width = 10, height = 9, standAlone = TRUE,
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
  par(mar = rep(0.2, 4))
  plot(NA, xlim = c(-4, n2+7), ylim = c(-5, n1+9), axes = FALSE,
       xlab = "", ylab = "")
  for (i in 1:n1){
    for (j in 1:n2){
      rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5, col = index[i,j], border = "white")
    }
  }
  lines(c(0.5, n2+0.5), c(0.5, 0.5))
  lines(c(0.5, 0.5), c(0.5, n1+0.5))
  nu_x=seq(20,80,by=10)
  for (i in 1:length(nu_x)){
    ind_nu = which.min(abs(nu - nu_x[i]))
    lines(c(0.5, -0.1), c(ind_nu, ind_nu))
    # text(-1+(i==1)*0.3, ind_nu+0.12, nu_x[i], cex = 1.5)
    text(-1, ind_nu+0.12, nu_x[i], cex = 1.5)
  }
  sig_x = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
  for (i in 1:length(sig_x)){
    ind_nu = which.min(abs(sig - sig_x[i]))
    lines(c(ind_nu, ind_nu), c(0.5, -0.1))
    text(ind_nu, -1.2, sig_x[i], cex = 1.5)
  }
  mtext("$\\sigma_1$", cex = 2, line = -2.3, side = 1)
  mtext("$\\nu_2$", cex = 2, line = -2.5, side = 2,las=2)
  all_cols = rev(c(col1, col2, col3))
  m = length(all_cols)
  vals=paste(round(rev(c(0,q1,q2,q3))*100, 2),"\\%", sep = "")
  placement = c(2, 4, 6, 8, 10)
  k = 35
  myseq = (100-k)/5 + seq(from = 0, to = k, length.out = m)
  rect(mean(myseq[3:4]), 64.2, mean(myseq[6:7]), 65+3.2, border = "grey70", col="grey70")
  text(mean(myseq[3:7]), 70, paste0("$99\\%$ Simulation error tolerance\nunder H$_0:\\alpha=0.05$"), col = "grey70", cex = 0.95)
  points(myseq, rep(63, m)+3.2, pch = 15,
         col = all_cols, cex = 6)
  for (i in 1:12){
    if (i > 3 && i < 7){
      text(myseq[i], 66.2, vals[i], cex = 0.78, col = "grey70")
    }else{
      text(myseq[i], 66.2, vals[i], cex = 0.78, col = "black")
    }
    
  }
  if (tikzPlot){
    dev.off()
  } else {
    print(paste0("Figure A", d_typ+7))
    # invisible(readline(prompt="Press [enter] to continue"))
  }
}
####################
# Figure A12 
####################
tost=pow_TOST_tiles[nu_range,,2]
atost=pow_aTOST_tiles[nu_range,,2]
xtostno=pow_xTOSTno_tiles[nu_range,,2]
xtostoff=pow_xTOSToff_tiles[nu_range,,2]
alpha = 0.05
tosth = hist(tost, breaks = 100, plot = F)
atosth = hist(atost, breaks = 100, plot = F)
xtostnoh = hist(xtostno, breaks = 100, plot = F)
xtostoffh = hist(xtostoff, breaks = 100, plot = F)
range_break = range(atosth$breaks,tosth$breaks,xtostoffh$breaks,xtostnoh$breaks)
breaks_=seq(range_break[1],range_break[2],l=70)
tosth = hist(tost, breaks = breaks_, plot = F)
atosth = hist(atost, breaks = breaks_, plot = F)
xtostoffh = hist(xtostoff, breaks = breaks_, plot = F)
xtostnoh = hist(xtostno, breaks = breaks_, plot = F)
lower = qbinom(.005,100000,alpha)/100000
upper = qbinom(.995,100000,alpha)/100000
gg_color_hue = function(n, alpha = 1) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100, alpha = alpha)[1:n]
}
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}
cols = gg_color_hue(5, alpha = 0.4)
cols = cols[-(1:2)]
cols = cols[c(3,2,1)]
cols[3] = t_col("chartreuse4", percent = 60)
mygrey = "grey50"
grey_transp = adjustcolor("grey30", alpha.f = 0.15)
cols = gg_color_hue(4)
if (tikzPlot){
  library(tikzDevice)
  tikz(paste0(main_dir, "A12_hist_size.tex"), width = 6.5, height = 5, standAlone = TRUE,
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
par(mar = c(5, 0, 2, 0))
plot(NA, xlim = c(-0.8, 6), ylim = c(0.9, 10.5), axes = FALSE, xlab = "", ylab = "")
max_height = 0.5*max(c(max(tosth$density), max(xtostoffh$density), max(atosth$density),  max(xtostnoh$density)))
K = length(tosth$density)
for (i in 1:K){
  rect(100*tosth$breaks[i], 8.5, 100*tosth$breaks[i+1], 8.5 + tosth$density[i]/max_height, col = cols[1], border = mygrey)
  rect(100*atosth$breaks[i], 6, 100*atosth$breaks[i+1], 6 + atosth$density[i]/max_height, col = cols[2], border = mygrey)
  rect(100*xtostnoh$breaks[i], 3.5, 100*xtostnoh$breaks[i+1], 3.5 + xtostnoh$density[i]/max_height, col = cols[3], border = mygrey)
  rect(100*xtostoffh$breaks[i], 1, 100*xtostoffh$breaks[i+1], 1 + xtostoffh$density[i]/max_height, col = cols[4], border = mygrey)
}
lines(c(-0.15,6), c(1,1))
lines(c(-0.15,6), c(3.5,3.5))
lines(c(-0.15,6), c(6,6))
lines(c(-0.15,6), c(8.5,8.5))
lines(c(5,5), c(-1, 11), lwd = 1.5, lty = 2)
rect(100*lower, -1, 100*upper, 11, col = grey_transp, border = NA)
axis(1, at = 0:6, cex.axis = 1.2)
axis(3, at = 0:6, cex.axis = 1.2)
mtext("Empirical Size (\\%)", side = 1, line = 2.65, adj = 0.6, cex = 1.25)
text(-0.15, 1, "cTOST*", pos = 2, cex = 1.15)
text(-0.15, 3.5, "cTOST", pos = 2, cex = 1.15)
text(-0.15, 6, "$\\alpha$-TOST", pos = 2, cex = 1.15)
text(-0.15, 8.5, "TOST", pos = 2, cex = 1.15)
text(3.75, 10, "99\\% Simulation Error\ntolerance under H$_0$: $\\alpha = 5\\%$", col = mygrey)
text(1.4, 10, "Approximately 33.3\\% of\nthe values are equal to 0", col = mygrey)
arrows(0.4, 10, 0.1, 10, length = 0.1, col = mygrey)
if (tikzPlot){
  dev.off()
} else {
  print(paste0("Figure A", 12))
  # invisible(readline(prompt="Press [enter] to continue"))
}
####################
# Figure A13
####################
tost=pow_TOST_tiles[nu_range,,1]
atost=pow_aTOST_tiles[nu_range,,1]
xtostno=pow_xTOSTno_tiles[nu_range,,1]
xtostoff=pow_xTOSToff_tiles[nu_range,,1]
tosth = hist(atost - tost, breaks = 30, plot = F)
xtostnoh = hist(atost - xtostno, breaks = 30, plot = F)
xtostoffh = hist(atost - xtostoff, breaks = 30, plot = F)
range_break = range(tosth$breaks,xtostnoh$breaks,xtostoffh$breaks)
breaks_=seq(range_break[1],range_break[2],l=50)
tosth = hist(atost - tost, breaks = breaks_, plot = F)
xtostnoh = hist(atost - xtostno, breaks = breaks_, plot = F)
xtostoffh = hist(atost - xtostoff, breaks = breaks_, plot = F)
gg_color_hue = function(n, alpha = 1) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100, alpha = alpha)[1:n]
}
t_col <- function(color, percent = 50, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  invisible(t.col)
}
cols = gg_color_hue(4)
cols = cols[-2]
mygrey = "grey50"
grey_transp = adjustcolor("grey30", alpha.f = 0.15)
if (tikzPlot){
  library(tikzDevice)
  tikz(paste0(main_dir, "A13_hist_power.tex"), width = 6.5, height = 5, standAlone = TRUE,
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
par(mar = c(4.3, 0, 2, 0))
plot(NA, xlim = c(-10, 18), ylim = c(0.9, 4.2), axes = FALSE, xlab = "", ylab = "")
max_height = 1.05*max(c(max(tosth$density), max(xtostnoh$density), max(xtostoffh$density)))
K = length(tosth$density)
for (i in 1:K){
  rect(100*tosth$breaks[i], 3, 100*tosth$breaks[i+1], 3 + tosth$density[i]/max_height, col = cols[1], border = mygrey)
  rect(100*xtostnoh$breaks[i], 2, 100*xtostnoh$breaks[i+1], 2 + xtostnoh$density[i]/max_height, col = cols[2], border = mygrey)
  rect(100*xtostoffh$breaks[i], 1, 100*xtostoffh$breaks[i+1], 1 + xtostoffh$density[i]/max_height, col = cols[3], border = mygrey)
}
axis(1, at = c(-4, 0, 4, 8, 12, 16), cex.axis = 1.2)
axis(3, at = c(-4, 0, 4, 8, 12, 16), cex.axis = 1.2)
mtext("Empirical Power Difference (\\%)", side = 1, line = 2.65, adj = 0.62, cex = 1.25)
text(-4.5, 3, "$\\alpha$-TOST$-$TOST", pos = 2, cex = 1.15)
text(-4.5, 2, "$\\alpha$-TOST$-$cTOST", pos = 2, cex = 1.15)
text(-4.5, 1, "$\\alpha$-TOST$-$cTOST*", pos = 2, cex = 1.15)
arrows(4, 4-0.15, 0.1, 4-0.15, length = 0.1, col = mygrey)
tt = round(mean((atost - tost)==0)*100, 1)
text(8.2, 4-0.15, paste0("Approximately ", tt, "\\% of\nthe values are equal to 0"), col = mygrey)
arrows(4, 4-1.3, 0.1, 4-1.3, length = 0.1, col = mygrey)
tt = round(mean((atost - xtostno)==0)*100, 1)
text(8.2, 4-1.3, paste0("Approximately ", tt, "\\% of\nthe values are equal to 0"), col = mygrey)
arrows(4, 4-2.3, 0.1, 4-2.3, length = 0.1, col = mygrey)
tt = round(mean((atost - xtostoff)==0)*100, 1)
text(8.2, 4-2.3, paste0("Approximately ", tt, "\\% of\nthe values are equal to 0"), col = mygrey)
lines(c(0,0), c(0.7, 4.4), lwd = 1.5, lty = 2)
mean((atost - tost)==0)
mean((atost - xtostno)==0)
if (tikzPlot){
  dev.off()
} else {
  print(paste0("Figure A", 13))
  # invisible(readline(prompt="Press [enter] to continue"))
}