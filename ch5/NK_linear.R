#======================================================
# Chapter 5, New Keynesian model (linear approximation)
#   modified on 2019/10/19
#======================================================

"%+%" <- function(x, y) paste(x, y, sep = "")

library(nleqslv)

#----------------------------------------

beta <- 0.99
mu <- 1.0
gamma <- 5
varrho <- 0.9
phi_pi <- 1.5
phi_y <- 0.5
rho_A <- 0.9
rho_v <- 0.7

kappa <- (1-varrho)*(1-varrho*beta)*(gamma+1)/varrho

para <- list(beta = beta, phi_pi = phi_pi, phi_y = phi_y, rho_A = rho_A, rho_v = rho_v, kappa = kappa)

#----------------------------------------

objfun0 <- function(X, maxT, X0,p = para) {
  
  nvar <- 2
  XX  <- matrix(X, maxT, nvar)
  x   <- c(XX[, 1], 0)
  ppi <- c(XX[, 2], 0)
  
  a  <- c(X0$a)
  v  <- c(X0$v)
  for(t in 1:maxT) {
    a[t+1] <- p$rho_A*a[t]
    v[t+1] <- p$rho_v*v[t]
  }
  ii <- p$phi_pi*ppi + p$phi_y*x + v
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- -x[t] + x[t+1]-(ii[t]-ppi[t+1])+(p$rho_A-1)*a[t]
    ret[t, 2] <- -ppi[t] + p$beta*ppi[t+1]+kappa*x[t]
  }
  return(c(ret))
}

rsltfun <- function(rslt, maxT, X0, Xss, p = para) {
  
  nvar <- 2
  XX  <- matrix(rslt, maxT, nvar)
  x   <- c(XX[, 1], 0)
  ppi <- c(XX[, 2], 0)
  
  a  <- c(X0$a)
  v  <- c(X0$v)

  for(t in 1:maxT) {
    a[t+1] <- p$rho_A*a[t]
    v[t+1] <- p$rho_v*v[t]
  }
  ii <- p$phi_pi*ppi + p$phi_y*x + v
  r <- ii-c(ppi[-1],0)
  
  dev <- list(x = x,  ppi = ppi, a = a, ii = ii, v = v, r = r)

  return(dev)
}

#----------------------------------------


maxT  <- 150
Xinit <- c(rep(0, maxT), rep(0, maxT))

a0a <- 0.01
v0a <- 0
X0a <- list(a = a0a, v = v0a)

objfun1 <- function(X) objfun0(X, maxT, X0a)
rslt1 <- nleqslv(Xinit, objfun1, control=list(maxit=10^4))
print(rslt1)

X1 <- rsltfun(rslt1$x, maxT, X0a)
X1[[1]]
X1[[2]]

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

maxT2 <- 50
plot(cbind(0:maxT2, c(0, X1$x[1:maxT2])), main = "x", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1$ii[1:maxT2])), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1$ppi[1:maxT2])), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1$r[1:maxT2])), main = "r", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1$a[1:maxT2])), main = "a", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1$v[1:maxT2])), main = "v", xlab = "", ylab = "", typ = "l")

#----------------------------------------

a0b <- 0
v0b <- 0.01
X0b <- list(a = a0b, v = v0b)

objfun2 <- function(X) objfun0(X, maxT, X0b)
rslt2 <- nleqslv(Xinit, objfun2, control=list(maxit=10^4))
print(rslt2)

X2 <- rsltfun(rslt2$x, maxT, X0b, ss)
X2[[1]]
X2[[2]]

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

maxT2 <- 50
plot(cbind(0:maxT2, c(0, X2$x[1:maxT2])), main = "x", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X2$ii[1:maxT2])), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X2$ppi[1:maxT2])), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X2$r[1:maxT2])), main = "r", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X2$a[1:maxT2])), main = "a", xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X2$v[1:maxT2])), main = "v", xlab = "", ylab = "", typ = "l")

