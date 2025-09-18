# rm(list=ls())
# B_aTOST = 10^5
require(tikzDevice)
require(ggplot2)
library(scales)
path = "./case_study/"
plot_path = paste0(path, "/tex")
options(scipen=999)
# construct ci
get_ci = function(dat,alpha=0.05){
  n = ncol(dat)
  m = apply(dat,1,mean)
  s = cov(t(dat))/n
  tval=qt(1-alpha,n-1)
  lower = m - sqrt(diag(s))*tval
  upper = m + sqrt(diag(s))*tval
  return(cbind(lower,upper))
}
# plot data
plot_dat_edited = function(dat,cte=log(1.25),xlim=NULL,main=NULL,var_nam=1:nrow(dat)){
  if(is.null(xlim)) xlim=range(dat) else xlim=xlim
  mylwd = 3
  p=nrow(dat)
  m = apply(dat,1,mean)
  ci = get_ci(dat)
  plot(rep(NA,p),p:1,xlim=xlim,yaxt='n',ylab="",xlab="",main=main)
  mtext("ECZ deposition", side=1, line = 2) # cex = cexsize,
  abline(h=1:p, v=seq(xlim[1], xlim[2], len=9), col="grey90", lwd=0.25, lty=2) # grid
  # axis(2,at=1:p,label=p:1)
  axis(2, at = 1:p, labels = rev(var_nam), las=2)
  lines(m,p:1,col="purple",lwd=mylwd)
  col_obs = alpha("red3",0.15)
  # for(i in 1:ncol(dat)) lines(dat[,i],p:1,col=col_obs,lwd=1.2)
  lines(ci[,2],p:1,lwd=mylwd)
  lines(ci[,1],p:1,lwd=mylwd)
  abline(v=c(-cte,cte),lwd=mylwd,lty=3,col="red3")
  legend("bottomleft",legend=c("Mean","90\\% CI", "Equivalence \n margins"),bty="n",
         lty=c(1,1,3),lwd=mylwd+1,pch=c(NA,NA,NA,NA),
         col=c("purple",1,"red3"),
         cex=0.85)
}
# load data
load(paste0(path, "porcine_skin.rda"))
dat = t(porcine_clean)
var_nam = c("\\textit{Stratum corneum} \n(0-20 $\\mu m$)",
            "\\textit{Viable epidermis} \n(20-160 $\\mu m$)",
            "\\textit{Upper dermis} \n(160-400 $\\mu m$)",
            "\\textit{Lower dermis} \n(400-800 $\\mu m$)")
if (tikzPlot){
  tikz(file=paste0(plot_path, "/app_raw_data.tex"),standAlone=TRUE,width=4.5,height=3.2)
}  
par(mfrow = c(1, 1),
    mai = c(0.1, 0.1, 0.1, 0.1),
    omi = c(0.5, 1.2, 0, 0))  # bltr
plot_dat_edited(dat, xlim=c(-1,1), var_nam=var_nam)
if (tikzPlot){
  dev.off()
}
#####################
# Analyze the data
#####################
alpha = 0.05
delta=log(1.25)
n=ncol(dat)
nu=n-1
p=nrow(dat)
m=apply(dat,1,mean)
s=cov(t(dat))/n
cat("Std. errors: \n", round(sqrt(diag(s)), 3), "\n")
cat("Means: \n", round(m, 3), "\n")
# TOST
source("aux_fun/atost_mvt.R")
tost_ci = ci_k(m, s, nu, alpha)
tost_ci_clean = round(tost_ci, 3)
tost_be = all(all(tost_ci[,1] > -delta), all(tost_ci[,2] < delta))
cat("TOST declares BE?", tost_be, "\n")
tost_res = cbind.data.frame(Lower=tost_ci[,1],
                            Mean=m,
                            Upper=tost_ci[,2])
# alpha-TOST
atost_out = get_alpha_TOST_MC_mv(B=B_aTOST,
                                     alpha=alpha,
                                     Sigma=s,
                                     nu=nu,
                                     delta=delta)
atost_alpha_star = atost_out$min
atost_ci = get_ci(dat,atost_alpha_star)
atost_ci_clean = round(atost_ci, 3)
atost_be = max(abs(atost_ci))<delta
cat("alpha-TOST declares BE?", atost_be, "\n")
atost_res = cbind.data.frame(Lower=atost_ci[,1],
                             Mean=m,
                             Upper=atost_ci[,2])
# cTOST
source("aux_fun/ctost_univ.R")
source("aux_fun/ctost_mvt.R")
ctost_out = get_ctost_mvt(alpha=alpha, Sigma=s, delta=rep(delta,p))
ctost_be = all(abs(m)<ctost_out$c_of_0)
cat("cTOST declares BE?", ctost_be, "\n")
ctost_ci_half_length = delta - ctost_out$c_of_0
ctost_ci = cbind(ci_L = m - ctost_ci_half_length, ci_U = ctost_ci_half_length + m)
# print(ctost_ci)
ctost_ci_clean = round(ctost_ci, 3)
ctost_res = cbind.data.frame(Lower=ctost_ci[,1],
                             Mean=m,
                             Upper=ctost_ci[,2])
# combine solutions
tost_row = paste0("[", tost_ci_clean[,1], ", ", tost_ci_clean[,2], "]")
atost_row = paste0("[", atost_ci_clean[,1], ", ", atost_ci_clean[,2], "]")
ctost_row = paste0("[", ctost_ci_clean[,1], ", ", ctost_ci_clean[,2], "]")
ci_df = rbind.data.frame(tost_row, atost_row, ctost_row)
rownames(ci_df) = c("Multivariate TOST", "Multivariate alpha-TOST ", "Multivariate cTOST")
colnames(ci_df) =  var_nam
ci_df
ci_df$"Bioequivalent?" = c(tost_be, atost_be, ctost_be)
ci_df
require(xtable)
if (tikzPlot){
  xtable(ci_df, type = "latex")
} else {
  print(xtable(ci_df, type = "latex"))
  print("Table 1")
  # invisible(readline(prompt="Press [enter] to continue"))
}
###################################################
# Figure 3: Confidence intervals & rejection region
###################################################
est_names = c("multivariate TOST", "multivariate $\\alpha$-TOST", "multivariate cTOST")
require(dplyr)
res_joined = full_join(full_join(cbind.data.frame(tost_res, id=rep(est_names[1], 4), x=rownames(tost_res)),
                                 cbind.data.frame(atost_res, id=rep(est_names[2], 4), x=rownames(atost_res))),
                       cbind.data.frame(ctost_res, id=rep(est_names[3], 4), x=rownames(ctost_res)))
res_joined$x = factor(res_joined$x, levels=rownames(tost_res))
res_joined$id = factor(res_joined$id, levels=est_names)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
color = gg_color_hue(4)[1:3]
# custom key to combine lines and points
draw_key_custom <- function(data, params, size, lwd=10) {
  # horizontal line with end caps
  line <- segmentsGrob(
    x0 = 0.1, x1 = 0.9,
    y0 = 0.5, y1 = 0.5,
    gp = gpar(col = data$colour, lwd = lwd)
  )
  # end caps (vertical lines)
  left_cap <- segmentsGrob(
    x0 = 0.1, x1 = 0.1,
    y0 = 0.4, y1 = 0.6,
    gp = gpar(col = data$colour, lwd = lwd)
  )
  right_cap <- segmentsGrob(
    x0 = 0.9, x1 = 0.9,
    y0 = 0.4, y1 = 0.6,
    gp = gpar(col = data$colour, lwd = lwd)
  )
  # center point
  point <- pointsGrob(
    x = 0.5, y = 0.5,
    pch = 16,
    gp = gpar(col = data$colour, size = data$size)
  )
  # all elements
  gTree(children = gList(line, left_cap, right_cap, point))
}
res_joined$id = factor(res_joined$id, levels = rev(c("multivariate TOST", "multivariate $\\alpha$-TOST", "multivariate cTOST")))
if (tikzPlot){
  tikz(file=paste0(plot_path, "/app_ci.tex"),standAlone=TRUE,width=13,height=6.5)
}
library(grid)
pp = ggplot(res_joined, aes(x, Mean)) +
  geom_rect(aes(xmin = 0.4, xmax = 4.6, ymin = -delta, ymax = delta),
            fill = alpha("gray", 0.03), color = NA) +
  geom_hline(yintercept=0, linetype="solid",
             color = "gray", size=2) +
  geom_point(aes(color = id), position = position_dodge(0.52), size=6,
             show.legend = T) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper, width=0.5, color=id),
                size=5, position = position_dodge(0.52),
                key_glyph = draw_key_custom,
                show.legend = T) +
  scale_y_continuous(limit = c(-0.55, 0.75), 
                     breaks = seq(-0.5, 0.75, len=6), 
                     labels = seq(-0.5, 0.75, len=6)) +
  coord_flip() +
  geom_hline(yintercept=-delta, linetype="dashed",
             color = "black", size=2) +
  geom_hline(yintercept=delta, linetype="dashed",
             color = "black", size=2) +
  theme_bw() +
  labs(y = "$\\theta$", x = NULL) +
  theme(text = element_text(size=30)) +
  scale_x_discrete(limits=rev, labels=rev(var_nam)) +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black")) +
  scale_color_manual("id",values=rev(color)) +
  theme(legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = c(0.4,-0.27))  +
  guides(colour = guide_legend(ncol = 3)) +
  guides(
    color = guide_legend(
      reverse = TRUE,
      override.aes = list(
        shape = 16,
        linewidth = 16,
        linetype = 1,
        segments.length = 14.5
      )
    )
  ) +
  theme(legend.key.spacing.x = unit(30, "pt")) +
  theme(panel.border = element_rect(colour = "black", fill=NA)) +
  theme(plot.margin = margin(1,1,0.1,0.1, "cm")) 
pp = pp + theme(plot.margin = unit(c(2,2.3,4,1), "lines")) + 
  annotation_custom(textGrob(label = "$-c_0$", hjust = 1, gp = gpar(cex = 2.5)),
                    ymin = -delta+0.04, ymax = -delta+0.04, 
                    xmin=4.8, xmax=4.8) + 
  annotation_custom(textGrob(label = "$c_0$", hjust = 1, gp = gpar(cex = 2.5)),
                    ymin = delta+0.025, ymax = delta+0.025, 
                    xmin=4.8, xmax=4.8) 
pp
gt <- ggplot_gtable(ggplot_build(pp))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)
if (tikzPlot){
  dev.off()
} else {
  print("Figure 3")
  # invisible(readline(prompt="Press [enter] to continue"))
}
