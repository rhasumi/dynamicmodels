#=======================================
# ノート5：MIUモデル
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
csvpath <- "F:/work/note/macro/R/ch5/"

library(nleqslv)

#----------------------------------------

beta <- 0.99
mu <- 1.0
gamma <- 1.0
rho_A <- 0.9
b <- 1
rho_psi <- 0.9
psi_ss <- 0
para <- list(beta = beta, mu = mu, gamma = gamma, rho_A = rho_A, b = b, rho_psi = rho_psi, psi_ss = psi_ss)

#----------------------------------------

findss <- function(p, Astar = 1) {
  istar <- 1/p$beta - 1;
  Cstar <- (1/(p$gamma+1)/p$mu)^(1/(p$gamma+1));
  pistar <- 0
  #mstar <- p$b*Cstar
  mstar <- p$b/(1-p$beta)*Cstar
  psistar <- p$psi_ss
  list(C = Cstar, ii = istar, A = Astar, ppi = pistar, m = mstar, psi = psistar)
}

objfun0 <- function(X, maxT, X0, Xss, p = para) {
  
  nvar <- 4
  XX  <- matrix(X, maxT, nvar)
  C   <- c(XX[, 1], Xss$C)
  A   <- c(X0$A   , XX[, 2])
  ppi <- c(XX[, 3], Xss$ppi)
  m   <- c(XX[, 4], Xss$m)
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- 1/m[t] - 1/C[t] + 0.99/(1+ppi[t+1])/C[t+1]
    ret[t, 2] <- -1 + (p$gamma+1)*p$mu*(C[t]/A[t])^(p$gamma+1)
    ret[t, 3] <- -(1+ppi[t+1])*m[t+1] + m[t]
    ret[t, 4] <- -log(A[t+1]) + p$rho_A*log(A[t])

  }
  return(c(ret))
}

#----------------------------------------

ss <- findss(para)

maxT  <- 100
Xinit <- c(rep(ss$C,maxT), rep(ss$A,maxT), 
           rep(ss$ppi,maxT), rep(ss$m,maxT))

Ainit1 <- 1.05
ppiinit1 <- 0
psiinit1 <- psi_ss
X0set1 <- list(A = Ainit1, ppi = ppiinit1, psi = psiinit1, m = ss$m)

objfun1 <- function(X) objfun0(X, maxT, X0set1, Xss = ss)
rslt1 <- nleqslv(Xinit, objfun1)

rsltfun <- function(rslt, maxT, X0, Xss, p = para) {
  XX <- matrix(rslt$x, maxT, 6)
  
  C   <- (c(Xss$C,  XX[, 1])/Xss$C -1)*100
  ii  <- (c(Xss$ii, XX[, 2])-Xss$ii)*100
  A   <- (c(Xss$A, X0$A, XX[1:(maxT-1), 3])/Xss$A-1)*100
  #ppi <- (c(Xss$ppi, X0$ppi, XX[1:(maxT-1), 4])-Xss$ppi)*100
  ppi <- (c(Xss$ppi, XX[, 4])-Xss$ppi)*100
  m   <- (c(Xss$m,  XX[, 5])/Xss$m -1)*100
  psi <- (c(Xss$psi, X0$psi, XX[1:(maxT-1), 6])-Xss$psi)*100
  
  list(C = C, ii = ii, A = A, ppi = ppi, m = m, psi = psi)
}

X1 <- rsltfun(rslt1, maxT, X0set1, ss)

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X1$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$m[1:51]), main = "m", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$ppi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X1$psi[1:51]), main = "psi", xlab = "", ylab = "", typ = "l")

#write.csv(X1, csvpath %+% "miu_techs.csv")

Ainit2 <- 1.0
ppiinit2 <- 0.0
psiinit2 <- -0.01
X0set2 <- list(A = Ainit2, ppi = ppiinit2, psi = psiinit2, m = ss$m)

objfun2 <- function(X) objfun0(X, maxT, X0set2, Xss = ss)
rslt2 <- nleqslv(Xinit, objfun2)

X2 <- rsltfun(rslt2, maxT, X0set2, ss)

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:50, X2$C[1:51]), main = "C", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$ii[1:51]), main = "i", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$m[1:51]), main = "m", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$ppi[1:51]), main = "pi", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$A[1:51]), main = "A", xlab = "", ylab = "", typ = "l")
plot(cbind(0:50, X2$psi[1:51]), main = "psi", xlab = "", ylab = "", typ = "l")

#write.csv(X2, csvpath %+% "miu_mss.csv")

