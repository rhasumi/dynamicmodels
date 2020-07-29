#=======================================
# Chapter 4, RBC model
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

library(nleqslv)

#----------------------------------------

alpha <- 0.3
beta <- 0.99
delta <- 0.025
mu <- 1.0
gamma <- 1.0
rho <- 0.9

para <- list(alpha = alpha, beta = beta, delta = delta, mu = mu, gamma = gamma, rho = rho)

#----------------------------------------

findss <- function(p, Astar = 1) {
  
  rstar <- 1/p$beta + p$delta - 1
  K_L <- (rstar/p$alpha/Astar)^(1/(p$alpha-1))
  Y_L <- Astar*K_L^p$alpha
  C_L <- Y_L-p$delta*K_L
  
  wstar <- (1-p$alpha)*Astar*K_L^p$alpha
  Lstar <- (wstar/(p$gamma+1)/p$mu)^(1/(p$gamma+1))*C_L^(-1/(p$gamma+1))
  Kstar <- K_L*Lstar
  Ystar <- Y_L*Lstar
  Cstar <- C_L*Lstar
  
  list(C = Cstar, L = Lstar, K = Kstar, Y = Ystar, w = wstar, r = rstar, A = Astar)
}



objfun0 <- function(X, maxT, X0, Xss, p = para) {
  
  nvar <- 7
  XX <- matrix(X, maxT, nvar)
  C  <- c(XX[, 1], Xss$C)
  L  <- c(XX[, 2], Xss$L)
  K  <- c(X0$K   , XX[, 3])
  Y  <- c(XX[, 4], Xss$Y)
  w  <- c(XX[, 5], Xss$w)
  r  <- c(XX[, 6], Xss$r)
  A  <- c(X0$A   , XX[, 7])
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- -w[t]/C[t]+(p$gamma+1)*p$mu*L[t]^p$gamma
    ret[t, 2] <- -C[t+1]/C[t]+p$beta*(r[t+1]-p$delta+1)
    ret[t, 3] <- -Y[t] + A[t]*K[t]^p$alpha*L[t]^(1-p$alpha)
    ret[t, 4] <- -w[t] + (1-p$alpha)*A[t]*K[t]^p$alpha*L[t]^(-p$alpha)
    ret[t, 5] <- -r[t]+p$alpha*A[t]*K[t]^(p$alpha-1)*L[t]^(1-p$alpha)
    ret[t, 6] <- -K[t+1]+Y[t]+(1-p$delta)*K[t]-C[t]
    ret[t, 7] <- -log(A[t+1])+p$rho*log(A[t])
  }
  return(c(ret))
}

#----------------------------------------

ss <- findss(para)

maxT  <- 150
Xinit <- c(rep(ss$C,maxT), rep(ss$L,maxT), rep(ss$K,maxT), rep(ss$Y,maxT), rep(ss$w,maxT), rep(ss$r,maxT), rep(ss$A,maxT))

K0 <- ss$K
A0 <- 1.01
X0 <- list(K = K0, A = A0)

objfun1 <- function(X) objfun0(X, maxT, X0, Xss = ss)
rslt1 <- nleqslv(Xinit, objfun1)

rsltfun <- function(rslt, maxT, X0, Xss, p = para) {
  XX <- matrix(rslt$x, maxT, 7)
  C0 <- c(Xss$C, XX[, 1])
  # K0 <- c(Xss$K, X0$K, XX[1:(maxT-1), 3])
  Y0 <- c(Xss$Y, XX[, 4])
  I0 <- Y0 - C0
  Iss <- Xss$Y - Xss$C
  
  sss <- 1-Xss$C/Xss$Y
  s <- I0/Y0
  
  C  <- (c(Xss$C, XX[, 1])/Xss$C-1)*100
  L  <- (c(Xss$L, XX[, 2])/Xss$L-1)*100
  K  <- (c(Xss$K, X0$K, XX[1:(maxT-1), 3])/Xss$K-1)*100
  Y  <- (c(Xss$Y, XX[, 4])/Xss$Y-1)*100
  w  <- (c(Xss$w, XX[, 5])/Xss$w-1)*100
  r  <- (c(Xss$r, XX[, 6])-Xss$r)*100
  A  <- (c(Xss$A, X0$A, XX[1:(maxT-1), 7])/Xss$A-1)*100
  
  I <- (I0/Iss-1)*100
  list(C = C, L = L, K = K, Y = Y, w = w, r = r, A = A, I = I, s = s)
}

X1 <- rsltfun(rslt1, maxT, X0, ss)

if (.Platform$OS.type == "windows") windows(7, 10)
par(ps = 14)
par(mfrow = c(4, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:maxT, X1$A), xlim = c(0,50), typ = "l", main = expression(A[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$Y), xlim = c(0,50), typ = "l", main = expression(Y[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$C), xlim = c(0,50), typ = "l", main = expression(C[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$K), xlim = c(0,50), typ = "l", main = expression(K[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$L), xlim = c(0,50), typ = "l", main = expression(L[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$I), xlim = c(0,50), typ = "l", main = expression(I[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$w), xlim = c(0,50), typ = "l", main = expression(w[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

plot(cbind(0:maxT, X1$r), xlim = c(0,50), typ = "l", main = expression(r[t]), xlab = "", ylab = "")
abline(h = 0, lty = 2, col = grey(0.5))

dev.copy2eps(file= figpath %+% "rbc1.eps")

if (.Platform$OS.type == "windows") windows(5*0.6, 3.55*0.6)
par(ps = 10*15/14)
par(mai = c(0.85, 0.8, 0.68, 0.35)*0.5)
plot(cbind(0:maxT, X1$s), xlim = c(0,50), typ = "l", main = expression(s[t]), xlab = "", ylab = "")
abline(h = 1-ss$C/ss$Y, lty = 2, col = grey(0.5))

dev.copy2eps(file= figpath %+% "rbc2s.eps")

