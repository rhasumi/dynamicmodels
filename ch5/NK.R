#=======================================
# Chapter 5, New Keynesian model
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")

library(nleqslv)

#----------------------------------------

beta <- 0.99
mu <- 1.0
gamma <- 5
eta <- 10
varrho <- 0.9
phi_pi <- 1.5
phi_y <- 0.5
rho_A <- 0.9
rho_v <- 0.7

para <- list(beta = beta, mu = mu, gamma = gamma, varrho = varrho, eta = eta, phi_pi = phi_pi, phi_y = phi_y, rho_A = rho_A, rho_v = rho_v)

#----------------------------------------

findss <- function(p) {
  pistar <- 0
  pitilstar <- 0
  vstar <- 0
  istar <- (1-p$beta)/p$beta
  Dstar <- 1/(1-p$varrho*p$beta)
  Fstar <- (p$eta-1)/p$eta/(1-p$varrho*p$beta)
  phistar <- (p$eta-1)/p$eta
  Astar <- 1
  Cstar <- ((p$eta-1)/p$eta/p$mu/(p$gamma+1))^(1/(p$gamma+1))
  
  list(C = Cstar, phi = phistar, ppi = pistar, ppitil = pitilstar,
       F = Fstar, D = Dstar,  A = Astar, ii = istar, v = vstar)
}

objfun0 <- function(X, maxT, X0, Xss, p = para) {
  
  nvar <- 4
  XX  <- matrix(X, maxT, nvar)
  C   <- c(XX[, 1], Xss$C)
  ppi <- c(XX[, 2], Xss$ppi)
  F   <- c(XX[, 3], Xss$F)
  D   <- c(XX[, 4], Xss$D)
  
  A  <- c(X0$A)
  v  <- c(X0$v)
  for(t in 1:maxT) {
    A[t+1] <- exp(p$rho_A*log(A[t]))
    v[t+1] <- p$rho_v*v[t]
  }
  phi <- p$mu*(C/A)^(p$gamma+1)*(p$gamma+1)
  ppitil <- p$eta/(p$eta-1)*F/D*(1+ppi)-1
  ii <- Xss$ii + p$phi_pi*ppi + p$phi_y*log(C/Xss$C/A) + v
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- -(1+ppi[t+1])*C[t+1]/C[t] + p$beta*(1+ii[t])
    ret[t, 2] <- -F[t] + phi[t] + p$varrho*p$beta*(1+ppi[t+1])^p$eta*F[t+1]
    ret[t, 3] <- -D[t] + 1 + p$varrho*p$beta*(1+ppi[t+1])^(p$eta-1)*D[t+1]
    ret[t, 4] <- -(1+ppi[t])^(1-p$eta) + (1-p$varrho)*(1+ppitil[t])^(1-p$eta)+p$varrho
  }
  return(c(ret))
}

rsltfun <- function(rslt, maxT, X0, Xss, p = para) {
  
  nvar <- 4
  XX  <- matrix(rslt, maxT, nvar)
  C   <- c(XX[, 1], Xss$C)
  ppi <- c(XX[, 2], Xss$ppi)
  F   <- c(XX[, 3], Xss$F)
  D   <- c(XX[, 4], Xss$D)
  
  A  <- c(X0$A)
  v  <- c(X0$v)
  for(t in 1:maxT) {
    A[t+1] <- exp(p$rho_A*log(A[t]))
    v[t+1] <- p$rho_v*v[t]
  }
  phi <- p$mu*(C/A)^p$gamma*(p$gamma+1)/A*C
  ppitil <- p$eta/(p$eta-1)*F/D*(1+ppi)-1
  ii <- Xss$ii + p$phi_pi*ppi + p$phi_y*log(C/Xss$C/A) + v
  
  ans0 <- list(C = C, phi = phi, ppi = ppi, ppitil = ppitil,
               F = F, D = D, A = A, ii = ii, v = v)
  
  C1   <- (c(Xss$C, C)/Xss$C-1)*100
  phi1 <- (c(Xss$phi, phi)/Xss$phi-1)*100
  ppi1 <- (c(Xss$ppi, ppi)-Xss$ppi)*100
  ppitil1 <- (c(Xss$ppitil, ppitil)-Xss$ppitil)*100
  F1 <- (c(Xss$F, F)/Xss$F-1)*100
  D1 <- (c(Xss$D, D)/Xss$D-1)*100
  A1 <- (c(Xss$A, A)/Xss$A-1)*100
  ii1 <- (c(Xss$ii, ii)-Xss$ii)*100
  v1 <- (c(Xss$v, v) - Xss$v)*100
  r1 <- c(Xss$ii-Xss$ppi,ii1[-length(ii1)]-ppi1[-1])
  
  ans1 <- list(C = C1, phi = phi1, ppi = ppi1, ppitil = ppitil1, 
       F = F1, D = D1, A = A1, ii = ii1, v = v1, r = r1)
  
  return(list(lev = ans0, dev = ans1))
}

#----------------------------------------

ss <- findss(para)

maxT  <- 150
Xinit <- c(rep(ss$C,maxT), rep(ss$ppi,maxT), rep(ss$F,maxT), rep(ss$D,maxT))

A0a <- 1.01
v0a <- 0
X0a <- list(A = A0a, v = v0a)

objfun1 <- function(X) objfun0(X, maxT, X0a, Xss = ss)
rslt1 <- nleqslv(Xinit, objfun1, control=list(maxit=10^4))
print(rslt1)

X1 <- rsltfun(rslt1$x, maxT, X0a, ss)
X1[[1]]
X1[[2]]

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X1[[2]]$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1[[2]]$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1[[2]]$ppi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1[[2]]$r[1:51]), main = "r", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1[[2]]$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1[[2]]$v[1:51]), main = "v", xlab = "", ylab = "", typ = "l")

write.csv(X1$dev, "nkm_techs.csv")


#----------------------------------------

A0b <- 1
v0b <- 0.01
X0b <- list(A = A0b, v = v0b)

objfun2 <- function(X) objfun0(X, maxT, X0b, Xss = ss)
rslt2 <- nleqslv(Xinit, objfun2, control=list(maxit=10^4))
print(rslt2)

X2 <- rsltfun(rslt2$x, maxT, X0b, ss)
X2[[1]]
X2[[2]]

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X2[[2]]$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2[[2]]$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2[[2]]$ppi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2[[2]]$r[1:51]), main = "r", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2[[2]]$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2[[2]]$v[1:51]), main = "v", xlab = "", ylab = "", typ = "l")

write.csv(X2$dev, "nkm_mps.csv")

