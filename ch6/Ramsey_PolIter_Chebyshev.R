#==================================================
# Chapter 6, Ramsey model (Chebyshev approximation)
#   modified on 2019/10/19
#==================================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

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

nt <- 20
at <- 0.5
bt <- 3.0

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

if (.Platform$OS.type == "windows") windows(7, 7)
par(ps = 20)
par(mai = c(0.85, 1, 0.25, 0.1))
plot(cbind(node, cc), typ = "l", xlab = expression(K[t]), ylab = expression(C[t]))
dev.copy2eps(file= figpath %+% "polit_ram_ct.eps")

plot(cbind(node, v), typ = "l", xlab = expression(K[t]), ylab = expression(V(K[t])))
dev.copy2eps(file= figpath %+% "polit_ram_v.eps")

#--------------------
# シミュレーション
#--------------------

K01 <- Kss*0.5
KK <- K01
CC <- c()

for(i in 1:31) {
  CC[i] <- hc(KK[i])
  KK[i+1] <- g(KK[i], CC[i])
}

par(ps = 15)
par(mai = c(0.85*0.75, 0.68*0.75, 0.68*0.75, 0.5))
par(mfrow = c(2,1))
plot(cbind(1:31, CC[1:31]), typ = "l", main = expression(C[t]), xlab = "", ylab = "")
axis(side=4, labels = c(expression(C[1]),expression("C*")), at = c(CC[1],Css), las = 2)

plot(cbind(1:31, KK[1:31]), typ = "l", main = expression(K[t]), xlab = "", ylab = "")
axis(side=4, labels = c(expression(K[1]),expression("K*")), at = c(KK[1],Kss), las = 2)

dev.copy2eps(file= figpath %+% "polit_ram_ck.eps")

