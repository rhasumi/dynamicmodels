#=======================================
# Chapter 5, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

fx <- function(x) sin(x)+1.5

windows(7, 7)
par(ps = 18)
par(mai = c(0.85, 0.85, 0.35, 0.35))
curve(fx, -0.5, 4, ylim = c(0, 3), ylab = "y")
abline(h = 0)
lines(c(0.5, 0.5), c(0, fx(0.5)))
lines(c(3.0, 3.0), c(0, fx(3.0)))
text(0.0, 0.75, "x = a")
text(3.5, 0.75, "x = b")
text(2.0, 2.7, "y = f(x)")
grids <- seq(0.5, 3.0, length =100)
dots <- rbind(c(0.5,0),cbind(grids, fx(grids)),c(3.0,0))
polygon(dots,col="gray")

dev.copy2eps(file= figpath %+% "integ.eps")

#---------------------------------------

nkm_techs <- read.csv("nkm_techs.csv")
cia_techs <- read.csv("cia_techs.csv")

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35*2)*0.5)

plot(cbind(0:50, nkm_techs$C[1:51]), col = 4, main = expression(C[t]), xlab = "", ylab = "", typ = "l", ylim = c(0,1))
lines(cbind(0:50, cia_techs$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ii[1:51]), col = 4, main = expression(i[t]), xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, cia_techs$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ppi[1:51]), col = 4, main = expression(pi[t]), xlab = "", ylab = "", typ = "l", ylim = c(-0.5, 0.1)/5)
axis(4, at = seq(-1, 1, 0.02), labels = seq(-1, 1, 0.02)*10)
lines(cbind(0:50, cia_techs$pi[1:51]/10), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$r[1:51]), col = 4, main = expression(r[t]), xlab = "", ylab = "", typ = "l", ylim = c(-0.5, 0.5)/5)
axis(4, at = seq(-1, 1, 0.05), labels = seq(-1, 1, 0.05)*10)
lines(cbind(0:50, cia_techs$r[1:51]/10), col = 2, lty = 2)

plot(cbind(0:50, cia_techs$m[1:51]), col = 2, lty = 2, main = expression(m[t]), xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, nkm_techs$A[1:51]), col = 4, main = expression(A[t]), xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, cia_techs$A[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_cia_techs.eps")

#---------------------------------------

nkm_mps <- read.csv("nkm_mps.csv")
cia_mss <- read.csv("cia_mss.csv")

#par(ps = 15)
#par(mfrow = c(3, 2))
#par(mai = c(0.85, 0.68, 0.68, 0.35*2)*0.5)

plot(cbind(0:50, nkm_mps$C[1:51]), col = 4, main = expression(C[t]), xlab = "", ylab = "", typ = "l", ylim = c(-1, 0.25))
lines(cbind(0:50, cia_mss$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ii[1:51]), col = 4, main = expression(i[t]), xlab = "", ylab = "", typ = "l", ylim = c(-4, 1)/5)
lines(cbind(0:50, cia_mss$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ppi[1:51]), col = 4, main = expression(pi[t]), xlab = "", ylab = "", typ = "l", ylim = c(-5, 0)/5)
lines(cbind(0:50, cia_mss$pi[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$r[1:51]), col = 4, main = expression(r[t]), xlab = "", ylab = "", typ = "l", ylim = c(0, 0.3))
lines(cbind(0:50, cia_mss$r[1:51]), col = 2, lty = 2)

plot(cbind(0:50, cia_mss$m[1:51]), col = 2, lty = 2, main = expression(m[t]), xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, nkm_mps$v[1:51]), col = 4, main = expression(paste(v[t], ", ", zeta[t])), xlab = "", ylab = "", typ = "l", ylim = c(-5, 5)/5)
lines(cbind(0:50, cia_mss$zeta[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_cia_ms.eps")

