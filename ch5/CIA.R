#=======================================
# Chapter 5, CIA model
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")

library(nleqslv)

#----------------------------------------

beta <- 0.99
mu <- 1.0
gamma <- 5 #1.0
rho_A <- 0.9
#rho_zeta <- 0.9
rho_zeta <- 0.7
zeta_ss <- 0
para <- list(beta = beta, mu = mu, gamma = gamma, rho_A = rho_A, rho_zeta = rho_zeta, zeta_ss = zeta_ss)

#----------------------------------------

findss <- function(p, Astar = 1) {
  pistar <- p$zeta_ss
  istar <- (1-p$beta+p$zeta_ss)/p$beta
  Cstar <- (p$beta/p$mu/(p$gamma+1)/(1+p$zeta_ss))^(1/(p$gamma+1))
  mstar <- (1+p$zeta_ss)*Cstar
  lamstar <- p$beta/(1+p$zeta_ss)/Cstar
  xistar <- (1+p$zeta_ss-p$beta)/(1+p$zeta_ss)/Cstar
  list(C = Cstar, ii = istar, A = Astar, pi = pistar, m = mstar, zeta = p$zeta_ss, lam = lamstar, xi = xistar)
}

objfun0 <- function(X, maxT, X0, Xss, p = para) {
  
  nvar <- 8
  XX  <- matrix(X, maxT, nvar)
  C   <- c(XX[, 1], Xss$C)
  ii  <- c(XX[, 2], Xss$ii)
  A   <- c(X0$A   , XX[, 3])
  pi <- c(XX[, 4], Xss$pi)
  m   <- c(XX[, 5], Xss$m)
  zeta <- c(X0$zeta , XX[, 6])
  lam <- c(XX[, 7], Xss$lam)
  xi <- c(XX[, 8], Xss$xi)
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- 1/C[t]- lam[t] - xi[t]
    ret[t, 2] <- A[t]*lam[t] - p$mu*(C[t]/A[t])^p$gamma*(p$gamma+1)
    if (t > 1) {
      ret[t, 3] <- -C[t] + m[t-1]/(1+pi[t])
      ret[t, 4] <- -(1+pi[t])*m[t] + (1+zeta[t])*m[t-1]
    } else {
      ret[t, 3] <- -C[t] + X0$m/(1+pi[t])
      ret[t, 4] <- -(1+pi[t])*m[t] + (1+zeta[t])*X0$m
    }
    ret[t, 5] <- p$beta*lam[t+1] + p$beta*xi[t+1] - (1+pi[t+1])*lam[t]
    ret[t, 6] <- p$beta*(ii[t]+1)*lam[t+1] - (1+pi[t+1])*lam[t]
    ret[t, 7] <- -log(A[t+1]) + p$rho_A*log(A[t])
    ret[t, 8] <- -(zeta[t+1]-p$zeta_ss) + p$rho_zeta*(zeta[t]-p$zeta_ss)
  }
  return(c(ret))
}

#----------------------------------------

ss <- findss(para)

maxT  <- 150
Xinit <- c(rep(ss$C,maxT), rep(ss$ii,maxT), rep(ss$A,maxT), 
           rep(ss$pi,maxT), rep(ss$m,maxT), rep(ss$zeta,maxT), rep(ss$lam,maxT), rep(ss$xi,maxT))

Ainit1 <- 1.01
zetainit1 <- zeta_ss
X0set1 <- list(A = Ainit1, zeta = zetainit1, m = ss$m)

objfun1 <- function(X) objfun0(X, maxT, X0set1, Xss = ss)
rslt1 <- nleqslv(Xinit, objfun1)

rsltfun <- function(rslt, maxT, X0, Xss, p = para) {
  XX <- matrix(rslt$x, maxT, 6)
  
  C   <- (c(Xss$C,  XX[, 1])/Xss$C -1)*100
  ii  <- (c(Xss$ii, XX[, 2])-Xss$ii)*100
  A   <- (c(Xss$A, X0$A, XX[1:(maxT-1), 3])/Xss$A-1)*100
  pi <- (c(Xss$pi, XX[, 4])-Xss$pi)*100
  m   <- (c(Xss$m,  XX[, 5])/Xss$m -1)*100
  zeta <- (c(Xss$zeta, X0$zeta, XX[1:(maxT-1), 6])-Xss$zeta)*100
  r <- c(Xss$ii-Xss$pi,ii[-length(ii)]-pi[-1])
  
  list(C = C, ii = ii, A = A, pi = pi, m = m, zeta = zeta, r = r)
}

X1 <- rsltfun(rslt1, maxT, X0set1, ss)

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X1$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$m[1:51]), main = "m", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$pi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
#plot(cbind(0:50, X1$zeta[1:51]), main = "zeta", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$r[1:51]), main = "r", xlab = "", ylab = "", typ = "l")

write.csv(X1, "cia_techs.csv")

Ainit2 <- 1.0
zetainit2 <- -0.01
X0set2 <- list(A = Ainit2, zeta = zetainit2, m = ss$m)

objfun2 <- function(X) objfun0(X, maxT, X0set2, Xss = ss)
rslt2 <- nleqslv(Xinit, objfun2)

X2 <- rsltfun(rslt2, maxT, X0set2, ss)

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X2$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$m[1:51]), main = "m", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$pi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
#plot(cbind(0:50, X2$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$zeta[1:51]), main = "zeta", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$r[1:51]), main = "r", xlab = "", ylab = "", typ = "l")

write.csv(X2, "cia_mss.csv")

