# rm(list=ls())
require(ggplot2)
mypath_output = "simulations/other_plots/"

###################
# uncomments to re-run (will take a while..)
###################
# mypath_output = "simulations/other_plots/A2_lambda_vs_rho_data"
# source("aux_fun/ctost_mvt.R")
# source("aux_fun/ctost_univ.R")
# rhotot = c(0, 0.5, 0.9)
# alphatot = c(0.05)
# for (rho in rhotot) {
#   for (alpha in alphatot) {
#     seed = 123
#     set.seed(seed)
#     sigmas = c(0.1, 0.1)
#     Sigma = build_Sigma(sigmas=sigmas,rho)
#     nu = 20
#     cte = log(1.25)
#     B = 10^5
#     B_aTOST = 10^4
#     thetas2_=log(seq(1,1.3,by=0.005))
#     thetas1_ = log(seq(1.3,1,by=-0.005))
#     tmp = matrix(NA,length(thetas1_),length(thetas2_))
#     colnames(tmp) = round(thetas2_, 4)
#     rownames(tmp) = round(thetas1_, 4)
#     pow_TOST_ = pow_alphaTOST_ = tmp
#     for (theta2 in thetas2_) {
#       for (theta1 in thetas1_) {
#         theta = c(theta1,theta2)
#         # TOST
#         pow_TOST_[which(thetas1_==theta1),which(thetas2_==theta2)] = 
#           power_TOST_MC_mv(alpha = alpha, theta = theta, Sigma = Sigma,
#                            nu = nu, delta = cte, B=B)$power_mult
#         sup_TOST = find_sup_x(alpha,Sigma,delta=cte,seed=10^5)
#       }
#     }
#     namefile = paste0(mypath_output, "corr_", rho, "_alpha_",  alpha,"_grid.RData")
#     save.image(file=namefile)
#   }
# }


####################
# plot 1
####################
source("aux_fun/ctost_mvt.R")
source("aux_fun/ctost_univ.R")
alpha = 0.05
nu = 20
cte = log(1.25)
# B = 10^5
rhotot = seq(0, 0.9999, len=1000)
sup_TOST = matrix(NA, length(rhotot), 2)
for (rho in rhotot) {
  
  sigmas = c(0.1, 0.1)
  Sigma = build_Sigma(sigmas=sigmas,rho)
  
  sup_TOST[which(rhotot %in% rho),] = find_sup_x(alpha,Sigma,
                                                 delta=cte,seed=12345)
  
}
df <- data.frame(
  rho = rep(rhotot, 2),
  lambda = c(sup_TOST[,1], sup_TOST[,2]),
  type = factor(rep(c("lambda1", "lambda2"), each = length(rhotot)))
)

col1 = c("#0072B2", "#0072B2")
# Create the plot
p1 = ggplot(df, aes(x = rho, y = lambda, color = type)) +
  geom_line(data = subset(df, type == "lambda1"), size = 1, linetype = "solid") +
  geom_line(data = subset(df, type == "lambda2"), size = 1, linetype = "dashed") +
  scale_color_manual(
    values = c("lambda1" = col1[1], "lambda2" = col1[2]),  # Blue and orange - colorblind friendly
    labels = c(expression(lambda[1]), expression(lambda[2]))
  ) +
  geom_abline(slope = cte, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, cte*1.1)) +
  labs(
    x = expression(rho),
    y = expression(lambda[j])
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.96, 0.04),  # Positioned to bottom right
    legend.justification = c(1, 0),   # Anchored at bottom right
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA)
  )
p1

#################################
# Plot 2
#################################

library(gridBase)
library(grid)

plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <-plotViewport(c(1.8,1,0,1))

require(ggplot2)

mypath_input = "simulations/other_plots/A2_lambda_vs_rho_data/"

rhotot   = 0.5
alphatot = 0.05
xseq        = seq(0, 0.25, by=0.05)
xseq[6]     = log(1.25)
xseq[7]     = 0.25
xseq_lab = expression("0", "0.05", "0.10", "0.15", "0.2", c[0], "0.25")

library(RColorBrewer)
my_orange_all = brewer.pal(5, "Oranges")
my_orange     = my_orange_all[3]

require(ggplot2)

pl_all = list()
uu=1

for (rho in rhotot) {
  for (alpha in alphatot) {
    
    namefile = paste0(mypath_input, "corr_", rho, "_alpha_",  alpha,"_grid.RData")
    load(file=namefile)
    
    # fix log scales used before
    colnames(pow_alphaTOST_) = colnames(pow_TOST_) = seq(thetas2_[1], thetas2_[length(thetas2_)], len=length(thetas2_))
    rownames(pow_alphaTOST_) = rownames(pow_TOST_) = seq(thetas1_[1], thetas1_[length(thetas1_)], len=length(thetas1_))
    
    library(reshape)
    rng  = c(0, max(pow_TOST_))
    dats = melt(t(pow_TOST_)) # use the TOST
    
    sqint       = seq(0, round(max(pow_TOST_)/5, 2)*5, length=20)
    dats$value2 = findInterval(dats$value, sqint, left.open=T)
    for (ii in 1:length(sqint)) {
      dats$value2[dats$value2==ii] = sqint[ii]
    }
    dats$Prob = dats$value2
    
    require(viridisLite)
    colors  = viridis(50)
    colors  = c(colors[round(seq(1,30,len=10))], colors[31:50])
    breaks2 = round(seq(0.2, 0, length=7), 2)
    
    xsup = sup_TOST[1]
    ysup = sup_TOST[2]
    
    my_lab = c(paste0("~lambda[1]==", round(xsup, 2)),
               paste0("~lambda[2]==", round(xsup, 1)))
    
    pl = ggplot(data = dats) + 
      geom_tile(aes(x=X1,y=X2,fill = Prob)) +
      scale_fill_gradientn(colours = colors, limits=range(breaks2),
                           breaks=breaks2, na.value="black",
                           guide = "legend") +
      geom_segment(x = log(1.25), xend = log(1.25), y = 0, yend = log(1.25),
                   linetype="dashed", color = my_orange, size=1) +
      geom_segment(y = log(1.25), yend = log(1.25), x = 0, xend = log(1.25),
                   linetype="dashed", color = my_orange, size=1) +
      scale_x_continuous(expand=c(0,0),
                         breaks=xseq, labels=xseq_lab) +
      scale_y_continuous(expand=c(0,0),
                         breaks=xseq, labels=xseq_lab) +
      coord_cartesian(xlim = range(thetas2_)-0.0025,
                      ylim = range(thetas1_)-0.0025,
                      clip = 'off') +
      xlab(expression(theta[1])) +
      ylab(expression(theta[2])) +
      theme(text = element_text(size=20),
            legend.key.size = unit(1, 'cm')) + 
      theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
      geom_segment(x = 0, xend = max(thetas1_)/(ysup/xsup), y = 0, yend = max(thetas1_),
                   linetype="solid", color = my_orange, size=1.5) + 
      geom_point(x=xsup,
                 y=ysup,
                 size=7,
                 shape=21,
                 stroke = 0.1,
                 col= alpha(my_orange, 0.25),
                 fill= paste0(my_orange,10)) +
      annotate("label", x = xsup+0.03, y = ysup-0.025, 
               label = sprintf("atop(~lambda[1] %s %s, lambda[2]==c[0])", 
                               "%~~%", round(xsup, 2)),
               parse = TRUE, size = 6) +
      theme(legend.justification="right",
            plot.margin = margin(1.51,0.755,0.1,1.5, "cm")) + # trbl
      guides(fill=guide_legend(title="Prob."))
    lab_offset = 0.02
    lab_size = 9
    pl = pl + annotate(geom = "text", y = max(thetas1_)+lab_offset, x = median(thetas1_), 
                       label = paste("~rho==~", rho), 
                       parse=T, size=lab_size,
                       color = "black")
    # pl
    pl_all[[uu]] = pl
    uu = uu+1
  }
}
p2 = pl_all[[1]]
pl_all = list(p1, p2)

# common theme to apply to both plots
common_theme <- theme_minimal() +
  theme(
    text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 15),
    legend.key.size = unit(1.1, "cm"),
    legend.key.height = unit(0.9, "cm"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    legend.title = element_text(size = 14),
    panel.grid = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
  )
# blue line values closest to x=0.5
blue_data <- df[df$type == "lambda1",]
# sort the data by rho to ensure proper interpolation
blue_data <- blue_data[order(blue_data$rho),]
# indices for points below and above x=0.5
below_idx <- max(which(blue_data$rho <= 0.5))
above_idx <- min(which(blue_data$rho >= 0.5))
# linear interpolation to find exact y value at x=0.5
intersection_y <- blue_data$lambda[below_idx] + 
  (blue_data$lambda[above_idx] - blue_data$lambda[below_idx]) * 
  (0.5 - blue_data$rho[below_idx]) / 
  (blue_data$rho[above_idx] - blue_data$rho[below_idx])
# update p1 with exact axis ranges and simplified x-axis labels
col1 = c("#0072B2", "#0072B2")
p1 = ggplot(df, aes(x = rho, y = lambda, color = type)) +
  geom_line(data = subset(df, type == "lambda1"), size = 1, linetype = "solid") +
  geom_line(data = subset(df, type == "lambda2"), size = 1, linetype = "dashed") +
  scale_color_manual(
    values = c("lambda1" = col1[1], "lambda2" = col1[2]),
    labels = c(expression(lambda[1]), expression(lambda[2]))
  ) +
  geom_abline(slope = cte, intercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0.5) +
  geom_point(
    data = data.frame(x = 0.5, y = intersection_y),
    aes(x = x, y = y), 
    inherit.aes = FALSE,  # avoid inheriting the color mapping
    color = my_orange,
    size = 7,
    shape = 21,
    stroke = 0.1,
    fill = paste0(my_orange, "80")
  ) +
  geom_point(
    data = data.frame(x = 0.5, y = log(1.25)),
    aes(x = x, y = y), 
    inherit.aes = FALSE,  # avoid inheriting the color mapping
    color = my_orange,
    size = 7,
    shape = 21,
    stroke = 0.1,
    fill = paste0(my_orange, "80")
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
    labels = c("0", "0.2", "0.4", "0.6", "0.8", "1.0"),  # labels with "0" instead of "0.0"
    minor_breaks = seq(0, 1, 0.1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 0.25),
    breaks = c(0, 0.05, 0.10, 0.15, 0.20, cte, 0.25),
    labels = c("0", "0.05", "0.10", "0.15", "0.20", expression(c[0]), "0.25"),
    minor_breaks = seq(0, 0.25, 0.025),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(rho),
    y = expression(lambda[k])
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 0.25), clip = "on") +
  common_theme +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.position = c(0.96, 0.04),
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    panel.grid.minor = element_line(color = "gray90", size = 0.1),
    panel.grid.major = element_line(color = "gray90", size = 0.2),
    panel.border = element_rect(color = "black", fill = NA)
  )
# update p2 (pl_all[[2]]) to match p1's aesthetics
p2 <- pl_all[[2]] +
  common_theme +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    legend.key.size = unit(0.8, 'cm'),
    panel.border = element_blank()
  )
# plots side by side with consistent spacing and legend
library(ggpubr)
plarr = ggarrange(p1, p2,
                  ncol = 2, nrow = 1,
                  common.legend = FALSE,
                  # labels = c("A", "B"),
                  font.label = list(size = 16, face = "bold"),
                  widths = c(1, 1.2))
if (tikzPlot){
  tit = paste0(mypath_output, "/A2_lambda_vs_rho_data.pdf")
  ggsave(file = tit, plot = plarr, width = 13, height = 6, dpi = 300)
} else {
  print(plarr)
  print("Figure A2")
  # invisible(readline(prompt="Press [enter] to continue"))
}
