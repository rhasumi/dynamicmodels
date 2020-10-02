#=======================================
# Chapter 7, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

beta <- 0.99
gamma <- 5 # ここ変更
varrho <- 0.9

kappa <- (1-varrho)*(1-varrho*beta)*(gamma+1)/varrho

phi_y1 <- -1
phi_y2 <- 3

phi_pi1 <- -(1-beta)*phi_y1/kappa+1
phi_pi2 <- -(1-beta)*phi_y2/kappa+1

if (.Platform$OS.type == "windows") windows(7, 7)
par(ps = 15)
par(mai = c(0.68*1.5, 0.68*1.5, 0.34, 0.34))

plot(NA, xlim = c(0, 2), ylim = c(0, 2), main = "", xlab = expression(phi[y]), ylab = expression(phi[pi]))
polygon(c(phi_y1, phi_y2, phi_y2, phi_y1, phi_y1), c(phi_pi1, phi_pi2, -1, -1, phi_pi1), col = "grey", border = "grey")
box()

text(1, 0.4, "indeterminate")
text(1, 1.5, "determinate")

dev.copy2eps(file= figpath %+% "nkm_indet.eps")


