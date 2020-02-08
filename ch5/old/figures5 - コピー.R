#=======================================
# ÉmÅ[Ég5: MIU/NKM
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "F:/work/hosei/note/figs/"
wd <- "F:/work/hosei/note/R/ch5/"

fx <- function(x) sin(x)+1.5

par(ps = 15)
curve(fx, -0.5, 4, ylim = c(0, 3), ylab = "y")
abline(h = 0)
lines(c(0.5, 0.5), c(0, fx(0.5)))
lines(c(3.0, 3.0), c(0, fx(3.0)))
text(0.2, 0.75, "x = a")
text(3.3, 0.75, "x = b")
text(2.0, 2.6, "y = f(x)")
grids <- seq(0.5, 3.0, length =100)
dots <- rbind(c(0.5,0),cbind(grids, fx(grids)),c(3.0,0))
polygon(dots,col="gray")

dev.copy2eps(file= figpath %+% "integ.eps")

#---------------------------------------

nkm_techs <- read.csv(wd %+% "nkm_techs.csv")
miu_techs <- read.csv(wd %+% "miu_techs.csv")

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, nkm_techs$C[1:51]), col = 4, main = "C", xlab = "", ylab = "", typ = "l", ylim = c(0, 5))
lines(cbind(0:50, miu_techs$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ii[1:51]), col = 4, main = "i", xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_techs$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ppi[1:51]), col = 4, main = "pi", xlab = "", ylab = "", typ = "l", ylim = c(-0.5, 0.5))
lines(cbind(0:50, miu_techs$ppi[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$phi[1:51]), col = 4, main = "phi", xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, miu_techs$m[1:51]), col = 2, lty = 2, main = "m", xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, nkm_techs$A[1:51]), col = 4, main = "A", xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_techs$A[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_miu_techs.eps")

#---------------------------------------

nkm_mps <- read.csv(wd %+% "nkm_mps.csv")
miu_mss <- read.csv(wd %+% "miu_mss.csv")

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, nkm_mps$C[1:51]), col = 4, main = "C", xlab = "", ylab = "", typ = "l", ylim = c(-5, 1))
lines(cbind(0:50, miu_mss$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ii[1:51]), col = 4, main = "i", xlab = "", ylab = "", typ = "l", ylim = c(-5, 5))
lines(cbind(0:50, miu_mss$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ppi[1:51]), col = 4, main = "pi", xlab = "", ylab = "", typ = "l", ylim = c(-5, 5))
lines(cbind(0:50, miu_mss$ppi[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$phi[1:51]), col = 4, main = "phi", xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, miu_mss$m[1:51]), col = 2, lty = 2, main = "m", xlab = "", ylab = "", typ = "l", ylim = c(-0.1, 0.1))

plot(cbind(0:50, nkm_mps$v[1:51]), col = 4, main = "v, psi", xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_mss$psi[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_miu_ms.eps")

#---------------------------------------

par(ps = 12)
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
plot(cbind(c(0.5, 4.5), c(0.5, 4.5)), typ = "l", xlim = c(-0.5, 5), ylim = c(-0.5, 5), xlab = "", ylab = "",  axes = F)
arrows(-0.5, 0, 5, 0, length = 0.15)
arrows(0, -0.5, 0, 5, length = 0.15, lty = 2)
lines(cbind(c(0.5, 4.5), c(4.0, 1.0)))
lines(cbind(c(2.5, 2.5), c(0, 2.5)), lty = 3)
lines(cbind(c(0, 2.5), c(2.5, 2.5)), lty = 3)
text(5, -0.5, label = expression(hat(x)[t]), cex = 1.5)
text(-0.5, 5, label = expression(pi[t+1]), cex = 1.5)
text(2.5, -0.5, label = 0, cex = 1.5)
text(-0.5, 2.5, label = expression(pi[ss]), cex = 1.5)

text(4.5, 1.5, label = expression(NKPC), cex = 1.5)
text(4.5, 4.0, label = expression(NKIS), cex = 1.5)

dev.copy2eps(file= figpath %+% "lognkm.eps")


