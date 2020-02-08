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

plot(cbind(0:50, nkm_techs$C[1:51]), col = 4, main = expression(C[t]), xlab = "", ylab = "", typ = "l", ylim = c(0, 5))
lines(cbind(0:50, miu_techs$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ii[1:51]), col = 4, main = expression(i[t]), xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_techs$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$ppi[1:51]), col = 4, main = expression(pi[t]), xlab = "", ylab = "", typ = "l", ylim = c(-0.5, 0.5))
lines(cbind(0:50, miu_techs$ppi[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_techs$phi[1:51]), col = 4, main = expression(phi[t]), xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, miu_techs$m[1:51]), col = 2, lty = 2, main = expression(m[t]), xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, nkm_techs$A[1:51]), col = 4, main = expression(A[t]), xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_techs$A[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_miu_techs.eps")

#---------------------------------------

nkm_mps <- read.csv(wd %+% "nkm_mps.csv")
miu_mss <- read.csv(wd %+% "miu_mss.csv")

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, nkm_mps$C[1:51]), col = 4, main = expression(C[t]), xlab = "", ylab = "", typ = "l", ylim = c(-5, 0))
lines(cbind(0:50, miu_mss$C[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ii[1:51]), col = 4, main = expression(i[t]), xlab = "", ylab = "", typ = "l", ylim = c(-5, 0))
lines(cbind(0:50, miu_mss$ii[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$ppi[1:51]), col = 4, main = expression(pi[t]), xlab = "", ylab = "", typ = "l")
lines(cbind(0:50, miu_mss$ppi[1:51]), col = 2, lty = 2)

plot(cbind(0:50, nkm_mps$phi[1:51]), col = 4, main = expression(phi[t]), xlab = "", ylab = "", typ = "l")

plot(cbind(0:50, miu_mss$m[1:51]), col = 2, lty = 2, main = expression(m[t]), xlab = "", ylab = "", typ = "l", ylim = c(-0.1, 0.1))

plot(cbind(0:50, nkm_mps$v[1:51]), col = 4, main = expression(paste(v[t], ", ", psi[t])), xlab = "", ylab = "", typ = "l", ylim = c(-5, 5))
lines(cbind(0:50, miu_mss$psi[1:51]), col = 2, lty = 2)

dev.copy2eps(file= figpath %+% "nkm_miu_ms.eps")





