#=======================================
# Chapter 3, Ramsey model with tax
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

library(nleqslv)

#----------------------------------------

alpha <- 0.3
beta <- 0.99
delta <- 0.25
At <- 1.0

#----------------------------------------

#Kt <- function(tauk) ((1/beta+delta-1)/((1-tauk)*alpha*At))^(1/(alpha-1))
#FCt <- function(K, g) (1-g)*At*K^alpha - delta*K

Kt <- function(tauk) ((1/beta+delta-1)/((1-tauk)*alpha*At))^(1/(alpha-1))
FCt <- function(K, g) At*K^alpha - delta*K - g

#----------------------------------------

objfun0 <- function(X, Z, maxT, X0, Xss = Css) {
  
  XX   <- matrix(X, maxT, 2)
  C    <- c(XX[, 1], Xss)
  K    <- c(X0, XX[, 2])
  tauc <- Z[[1]]
  tauk <- Z[[2]]
  g    <- Z[[3]]
  
  ret <- matrix(NA, maxT, 2)
  for(t in 1:maxT) {
    ret[t, 1] <- (1+tauc[t+1])*C[t+1]/(1+tauc[t])/C[t]-beta*((1-tauk[t+1])*alpha*At*K[t+1]^(alpha-1)-delta+1)
    ret[t, 2] <- K[t+1]-At*K[t]^alpha-(1-delta)*K[t]+C[t]+g[t]
  }
  return(c(ret))
}

Kt0 <- Kt(0)
FCt0 <- FCt(Kt0, 0.1)

#----------------------------------------

#--------------------
# シミュレーション１
#--------------------

maxT  <- 30

tauc1 <- rep(0, maxT+1)
tauk1 <- rep(0, maxT+1)

for(t in 10:(maxT+1))
  tauc1[t] <- 0.1

g <- rep(0.1, maxT+1)

Z1 <- list(tauc1, tauk1, g)

Kss0 <- Kt(tauk = tauk1[1])
Css0 <- FCt(Kss0, g = g[1])
Xinit1 <- c(rep(Css0,maxT), rep(Kss0,maxT))

Kss1 <- Kt(tauk = tauk1[maxT+1])
Css1 <- FCt(Kss1, g = g[maxT+1])

objfun1 <- function(X) objfun0(X, Z1, maxT, Kss0, Css1)
rslt1 <- nleqslv(Xinit1, objfun1)
rslt1x <- matrix(rslt1$x, maxT, 2)

Ct1 <- c(Css0, rslt1x[, 1],  Css1)
Kt1 <- c(Kss0, Kss0, rslt1x[, 2])

par(ps = 15)
par(mai = c(0.85*0.75, 0.68*0.75, 0.68*0.75, 0.35*0.75))
par(mfrow = c(2,1))
plot(cbind(0:(maxT+1), Ct1), typ = "l", main = expression(C[t]), xlab = "", ylab = "")
segments(0, Css0, maxT+1, Css0, lty = 2, col = "blue")

plot(cbind(0:(maxT+1), Kt1), typ = "l", main = expression(K[t]), xlab = "", ylab = "")
segments(0, Kss0, maxT+1, Kss0, lty = 2, col = "blue")

dev.copy2eps(file= figpath %+% "pubck1.eps")

#--------------------
# シミュレーション２
#--------------------

tauc2 <- rep(0, maxT+1)
tauk2 <- rep(0, maxT+1)

for(t in 10:(maxT+1))
  tauk2[t] <- 0.1

Z2 <- list(tauc2, tauk2, g)
Xinit2 <- c(rep(Css0,maxT), rep(Kss0,maxT))

Kss2 <- Kt(tauk = tauk2[maxT+1])
Css2 <- FCt(Kss2, g = g[maxT+1])

objfun2 <- function(X) objfun0(X, Z2, maxT, Kss0, Css2)
rslt2 <- nleqslv(Xinit2, objfun2)
rslt2x <- matrix(rslt2$x, maxT, 2)

Ct2 <- c(Css0, rslt2x[, 1],  Css2)
Kt2 <- c(Kss0, Kss0, rslt2x[, 2])

par(ps = 15)
par(mai = c(0.85*0.75, 0.68*0.75, 0.68*0.75, 0.35*0.75))
par(mfrow = c(2,1))
plot(cbind(0:(maxT+1), Ct2), typ = "l", main = expression(C[t]), xlab = "", ylab = "")
segments(0, Css0, maxT+1, Css0, lty = 2, col = "blue")

plot(cbind(0:(maxT+1), Kt2), typ = "l", main = expression(K[t]), xlab = "", ylab = "")
segments(0, Kss0, maxT+1, Kss0, lty = 2, col = "blue")

dev.copy2eps(file= figpath %+% "pubck2.eps")

