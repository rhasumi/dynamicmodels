#=======================================
# Chapter 2, Figures of appendix
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

beta <- 0.99
alpha <- 0.3
delta <- 0.25

Kss <- ((1/beta+delta-1)/alpha)^(1/(alpha-1))
Css <- Kss^alpha - delta*Kss

U <- function (C) log(C)
g <- function (K, C) K^alpha + (1-delta)*K - C

y <- 0
for(i in 0:200)
  y <- y + beta^i * U(Css)

#---------------------------------------
# チェビシェフ近似補間（1次）
#---------------------------------------

library(nleqslv) # nleqslv

xj <- function(n = 10){
  -cos(0.5*(2*(1:n)-1)*pi/n)
}

trans0 <- function(x, a = -1, b = 1) {
  y <- (x-a)/(b-a)*2-1
  y
}

trans <- function(x, a = -1, b = 1) {
  y <- (x-a)/(b-a)*2-1
  y[which(y < -1)] <- -1
  y[which(y > 1)] <- 1
  y
}

rev <- function(y, a = -1, b = 1) (y+1)*(b-a)*0.5+a

fn <- function(x, cc, nt){
  size <- length(x)
  Tk <- matrix(1, size[1], 1)
  Tk <- cbind(Tk, x)
  for (k in 2:(nt-1))
    Tk <- cbind(Tk, 2*x*Tk[, k]-Tk[, k-1])
  Tk %*% cc
}

cheb <- function(x, coef, a, b, n) {
  fn(trans(matrix(x, 1), a, b), coef, n) 
}

#---------------------------------------
# Policy Iteration
#---------------------------------------

nt <- 40
at <- 0.05
bt <- 4.0

node <- rev(xj(nt), at[1], bt[1])
apxfun <- function(x, coef) fn(trans(x, at, bt), coef, nt)

niter <- 1000
eps <- 0.0001

# initialization
f <- function (x) -45 + sqrt(x)
cc <- rep(0, nt)
v <- rep(0, nt)
coef <- rep(0, nt)
coef2 <- rep(0, nt)

for(iter in 1:niter) {
  print(iter)
  v1 <- c(v)
  c1 <- c(cc)
  coef1 <- coef
  
  for(i in 1:nt) {
    obj <- function(ct) {
      st <- node[i]
      sa <- matrix(g(st,ct), 1)
      return(-U(ct)-beta*f(sa))
    }
    opt <- optimize(obj, c(0.001,  min((1 - delta) *node[i]+node[i]^alpha,bt*0.99)))
    cc[i] <- opt$minimum
  }
  fobj1 <- function(coef) cc - apxfun(node, coef)
  coef1 <- nleqslv(coef1, fobj1)$x
  hc <- function(x) cheb(x, coef1, at, bt, nt)
  
  for(i in 1:nt) {
    st <- node[i]
    ct <- max(min(hc(st), st), 0.001)
    y <- U(ct)
    for(k in 1:1000) {
       st <- g(st, ct)
       ct <- max(min(hc(st), st), 0.001)
       y <- y + beta^k * U(ct)
      }
    v[i] <- y
  }
  if (identical(abs(c(cc) - c1) < eps ,rep(T, nt))) break
  
  fobj2 <- function(coef) v - apxfun(node, coef)
  coef2 <- nleqslv(coef2, fobj2)$x
  f <- function(x) cheb(x, coef2, at, bt, nt)
  
  plot(cbind(node, v), typ = "l")
  print(abs(c(cc) - c1))
}

#---------------------------------------
# Plot
#---------------------------------------

#deltaK <- function(K) K^alpha - delta*K
#y1 <- c()
#for(i in 1:nt) 
#  y1[i] <- deltaK(node[i])
#
#windows(7, 5)
#par(mai = c(0.85, 0.88,0.85, 0.35))
#par(ps = 15)
#plot(cbind(node, cc), typ = "l", xlab = expression(K[t]), ylab = expression(C[t]), xlim = c(0,2.5), main = "Policy Function")
#points(cbind(Kss, Css), pch = 19)
#text(cbind(Kss+0.20, Css-0.1), expression(paste("(",K[ss],", ", C[ss],")")))
#lines(cbind(node, y1), col = 2, lty = 2)
#
#legend(cbind(1.8, 1.4), legend = c(expression(paste(C[t]," = F(", K[t], ")")), expression(paste(Delta,K[t]," = 0"))), col = 1:2, lty = 1:2)
#
#text(cbind(0.5, 1.2), expression(paste(Delta,K[t]," < 0")), col = 2)
#text(cbind(2.0, 0.4), expression(paste(Delta,K[t]," > 0")), col = 2)
#
# dev.copy2eps(file= figpath %+% "polit_ram_ct2.eps")


#--------------------
# シミュレーション
#--------------------

K01 <- Kss*0.5
Kt1 <- K01
Ct1 <- c()

for(i in 1:31) {
  Ct1[i] <- hc(Kt1[i])
  Kt1[i+1] <- g(Kt1[i], Ct1[i])
}

K02 <- Kss*2
Kt2 <- K02
Ct2 <- c()
for(i in 1:31) {
  Ct2[i] <- hc(Kt2[i])
  Kt2[i+1] <- g(Kt2[i], Ct2[i])
}

#--------------------
# 位相図の作図
#--------------------

At <- 1
Kt <- ((1/beta+delta-1)/alpha/At)^(1/(alpha-1))
FCt <- function(K) At*K^alpha - delta*K

maxT <- 31

windows(7, 7)
par(ps = 15)
par(mai = c(0.85, 0.88, 0.35, 0.35))
par(mfrow = c(1,1))
#curve(FCt, 0, 3, ylim = c(0.3, 1.7), xlim = c(0.05, 2.5), xlab = expression(K[t]), ylab = expression(C[t]), main = "Phase Diagram",axes = F)
curve(FCt, 0, 3, ylim = c(0.3, 1.3), xlim = c(0.05, 2.5), xlab = expression(K[t]), ylab = expression(C[t]), main = "",axes = F)
abline(h = 0)
lines(cbind(c(Kt, Kt),c(0, 100)))
points(Kt,FCt(Kt), pch=19)

axis(2)
axis(1, labels = c(0:3), at = 0:3)
abline(v = c(K01,K02), lty = 2)
box()

axis(1, labels = c(expression(K[1]^A)), at = K01)
axis(1, labels = c(expression(K[1]^B)), at = K02)
axis(1, labels = c(expression("K*")), at = Kss)

for(i in 1:12){
  arrows(Kt1[i], Ct1[i], Kt1[i+1], Ct1[i+1], length = 0.05)
  arrows(Kt2[i], Ct2[i], Kt2[i+1], Ct2[i+1], length = 0.05)
}
for(i in 13:(maxT-1)){
  segments(Kt1[i], Ct1[i], Kt1[i+1], Ct1[i+1])
  segments(Kt2[i], Ct2[i], Kt2[i+1], Ct2[i+1])
}


FCa <- function(K, C) beta*(alpha*g(K,C)^(alpha-1)-delta+1)*C

Kt1a <- K01
Ct1a <- hc(K01)*1.05
for(i in 1:31) {
  Kt1a[i+1] <- g(Kt1a[i], Ct1a[i])
  Ct1a[i+1] <- FCa(Kt1a[i], Ct1a[i])
}

Kt1b <- K01
Ct1b <- hc(K01)*0.95
for(i in 1:31) {
  Kt1b[i+1] <- g(Kt1b[i], Ct1b[i])
  Ct1b[i+1] <- FCa(Kt1b[i], Ct1b[i])
}

Kt2a <- K02
Ct2a <- hc(K02)*1.05
for(i in 1:31) {
  Kt2a[i+1] <- g(Kt2a[i], Ct2a[i])
  Ct2a[i+1] <- FCa(Kt2a[i], Ct2a[i])
}

Kt2b <- K02
Ct2b <- hc(K02)*0.95
for(i in 1:31) {
  Kt2b[i+1] <- g(Kt2b[i], Ct2b[i])
  Ct2b[i+1] <- FCa(Kt2b[i], Ct2b[i])
}

for(i in 1:maxT){
  if (!is.nan(Kt1a[i]))
    arrows(Kt1a[i], Ct1a[i], Kt1a[i+1], Ct1a[i+1], length = 0.05, col = 2, lty = 1)
  if (!is.nan(Kt1b[i]))
    arrows(Kt1b[i], Ct1b[i], Kt1b[i+1], Ct1b[i+1], length = 0.05, col = 4, lty = 1)
  if (!is.nan(Kt2a[i]))
    arrows(Kt2a[i], Ct2a[i], Kt2a[i+1], Ct2a[i+1], length = 0.05, col = 2, lty = 1)
  if (!is.nan(Kt2b[i]))
    arrows(Kt2b[i], Ct2b[i], Kt2b[i+1], Ct2b[i+1], length = 0.05, col = 4, lty = 1)
}

points(cbind(Kt1a[1], Ct1a[1]), pch = 15, col = 2)
points(cbind(Kt1[1], Ct1[1]), pch = 4)
points(cbind(Kt1b[1], Ct1b[1]), pch = 17, col = 4)

points(cbind(Kt2a[1], Ct2a[1]), pch = 15, col = 2)
points(cbind(Kt2[1], Ct2[1]), pch = 4)
points(cbind(Kt2b[1], Ct2b[1]), pch = 17, col = 4)


maxT2 <- 100

Kt3a <- Kss-0.01
Ct3a <- Css+0.01
for(i in 1:maxT2) {
  Kt3a[i+1] <- g(Kt3a[i], Ct3a[i])
  Ct3a[i+1] <- FCa(Kt3a[i], Ct3a[i])
}
end <- min(which(is.nan(Ct3a)))
#Kt3a[end] <- Kt1a[8]+0.04
#Ct3a[end] <- 1.6

Kt3b <- Kss+0.01
Ct3b <- Css-0.01
for(i in 1:maxT2) {
  Kt3b[i+1] <- g(Kt3b[i], Ct3b[i])
  Ct3b[i+1] <- FCa(Kt3b[i], Ct3b[i])
}

for(i in 1:maxT2){
  if (!is.nan(Kt3a[i]))
    arrows(Kt3a[i], Ct3a[i], Kt3a[i+1], Ct3a[i+1], length = 0.05, col = 1, lty = 1)
  if (!is.nan(Kt3b[i]))
    arrows(Kt3b[i], Ct3b[i], Kt3b[i+1], Ct3b[i+1], length = 0.05, col = 1, lty = 1)
}

text(Kss+0.07, Css+0.07, "F")

dev.copy2eps(file= figpath %+% "phase_ram_ck.eps")

