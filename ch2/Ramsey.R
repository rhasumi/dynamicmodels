#=======================================
# Chapter 2, Ramsey model
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../"

library(nleqslv)

#----------------------------------------

alpha <- 0.3
beta <- 0.99
delta <- 0.25
At <- 1.0

#----------------------------------------

Kt <- ((1/beta+delta-1)/alpha/At)^(1/(alpha-1))
FCt <- function(K) At*K^alpha - delta*K
Kmax <- (delta/At)^(1/(alpha-1))

windows(7, 7)
par(ps = 15)
par(mai = c(0.85, 0.85, 0.35, 0.35))
curve(FCt, 0, Kmax, ylim = c(0, 1), xlab = expression(K), ylab = expression(C))
abline(h = 0)

lines(cbind(c(Kt, Kt),c(0, 100)))
points(Kt,FCt(Kt), pch=19)

text(4, 0.4, expression(paste(Delta,K[t]," = 0",sep="")))
text(2, 0.95, expression(paste(Delta,C[t]," = 0",sep="")))

axis(1, labels = c(expression(K[max])), at = Kmax)
segments(Kmax, -1, Kmax, 0, lty = 2)

dev.copy2eps(file= figpath %+% "dcdk.eps")

#----------------------------------------

Kss <- ((1/beta+delta-1)/alpha/At)^(1/(alpha-1))
Css <- At*Kss^alpha - delta*Kss

objfun0 <- function(X, maxT, X0, Xss = Css) {
  
  XX  <- matrix(X, maxT, 2)
  C   <- c(XX[, 1], Xss)
  K   <- c(X0, XX[, 2])
  
  ret <- matrix(NA, maxT, 2)
  for(t in 1:maxT) {
    ret[t, 1] <- C[t+1]/C[t]-beta*(alpha*At*K[t+1]^(alpha-1)-delta+1)
    ret[t, 2] <- K[t+1]-At*K[t]^alpha-(1-delta)*K[t]+C[t]
  }
  return(c(ret))
}

maxT  <- 30
Xinit <- c(rep(Css,maxT), rep(Kss,maxT))

#--------------------
# シミュレーション１
#--------------------

K01 <- Kss*0.5
objfun1 <- function(X) objfun0(X, maxT, K01)
rslt1 <- nleqslv(Xinit, objfun1)

rslt1x <- matrix(rslt1$x, maxT, 2)

Ct1 <- c(rslt1x[, 1],  Css)
Kt1 <- c(K01, rslt1x[, 2])

#--------------------
# シミュレーション２
#--------------------

K02 <- Kss*2
objfun2 <- function(X) objfun0(X, maxT, K02)
rslt2 <- nleqslv(Xinit, objfun2)

rslt2x <- matrix(rslt2$x, maxT, 2)

Ct2 <- c(rslt2x[, 1],  Css)
Kt2 <- c(K02, rslt2x[, 2])

#--------------------
# 作図
#--------------------

par(ps = 15)
par(mai = c(0.85*0.75, 0.68*0.75, 0.68*0.75, 0.35*0.75))
par(mfrow = c(2,1))
plot(cbind(1:(maxT+1), Ct1), typ = "l", main = expression(C[t]), xlab = "", ylab = "")
plot(cbind(1:(maxT+1), Kt1), typ = "l", main = expression(K[t]), xlab = "", ylab = "")

dev.copy2eps(file= figpath %+% "ramck.eps")

par(ps = 15)
par(mai = c(0.85, 0.85, 0.35, 0.35))
par(mfrow = c(1,1))
curve(FCt, 0, 3, ylim = c(0.5, 1.1), xlab = expression(K[t]), ylab = expression(C[t]))
abline(h = 0)
lines(cbind(c(Kt, Kt),c(0, 100)))
points(Kt,FCt(Kt), pch=19)

for(i in 1:12){
  arrows(Kt1[i], Ct1[i], Kt1[i+1], Ct1[i+1], length = 0.1)
  arrows(Kt2[i], Ct2[i], Kt2[i+1], Ct2[i+1], length = 0.1)
}
for(i in 13:(maxT-1)){
  segments(Kt1[i], Ct1[i], Kt1[i+1], Ct1[i+1])
  segments(Kt2[i], Ct2[i], Kt2[i+1], Ct2[i+1])
}

dev.copy2eps(file= figpath %+% "dcdk2.eps")

