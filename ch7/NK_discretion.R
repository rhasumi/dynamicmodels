#=======================================
# Chapter 7, Optimal MP (discretion)
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

library(nleqslv)

#----------------------------------------

beta <- 0.99
gamma <- 5
eta <- 10
varrho <- 0.9
phi_pi <- 1.5
phi_y <- 0.5

kappa <- (1-varrho)*(1-varrho*beta)*(gamma+1)/varrho
lambda <- (1-varrho)*(1-varrho*beta)/eta/varrho*(1+gamma)

para <- list(beta = beta, kappa = kappa, lambda = lambda)

#----------------------------------------

objfun0 <- function(X, maxT, X0, p = para) {
  
  nvar <- 3
  XX  <- matrix(X, maxT, nvar)
  x   <- c(XX[, 1], 0)
  ppi <- c(XX[, 2], 0)
  ii <- c(XX[, 3], 0)
  
  e  <- c(X0$e, rep(0, maxT))
  
  ret <- matrix(NA, maxT, nvar)
  for(t in 1:maxT) {
    ret[t, 1] <- -ppi[t] + p$beta*ppi[t+1]+kappa*x[t]+e[t]
    ret[t, 2] <- -x[t] + x[t+1]-(ii[t]-ppi[t+1])
    ret[t, 3] <- p$lambda*x[t]+p$kappa*ppi[t]
  }
  return(c(ret))
}

rsltfun <- function(rslt, maxT, X0, p = para) {
  
  nvar <- 3
  XX  <- matrix(rslt, maxT, nvar)
  x   <- c(XX[, 1], 0)
  ppi <- c(XX[, 2], 0)
  ii <- c(XX[, 3], 0)
  
  e  <- c(X0$e, rep(0, maxT))
  r <- ii-c(ppi[-1],0)

  dev <- list(x = x,  ppi = ppi, ii = ii, e = e, r = r)
  
  cumbeta <- c(1, cumprod(rep(p$beta, maxT)))
  W <- -0.5*sum(cumbeta*ppi^2+p$lambda*cumbeta*x^2)


  return(list(dev = dev, W = W))
}

#----------------------------------------

maxT  <- 150
Xinit <- c(rep(0, maxT), rep(0, maxT), rep(0, maxT))

e0 <- 1
X0 <- list(e = e0)

objfun1 <- function(X) objfun0(X, maxT, X0)
rslt1 <- nleqslv(Xinit, objfun1, control=list(maxit=10^4))
print(rslt1)

X1 <- rsltfun(rslt1$x, maxT, X0)
X1[[1]]
X1[[2]]/10^4

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

maxT2 <- 10
plot(cbind(0:maxT2, c(0, X1[[1]]$x[1:maxT2])), main = expression(x[t]),
  xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1[[1]]$ii[1:maxT2])), main = expression(i[t]),
  xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1[[1]]$ppi[1:maxT2])), main = expression(pi[t]),
  xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, -X1[[1]]$ppi[1:maxT2])), main = expression(psi[t]),
  xlab = "", ylab = "", typ = "l")
plot(cbind(0:maxT2, c(0, X1[[1]]$e[1:maxT2])), main = expression(e[t]),
  xlab = "", ylab = "", typ = "l")

write.csv(X1[[1]], "nkm_discretion.csv", row.names = F)

Y1 <- read.csv("nkm_commit.csv")

par(ps = 15)
par(mfrow = c(3, 2))
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

maxT2 <- 10
plot(cbind(0:maxT2, c(0, X1[[1]]$x[1:maxT2])), main = expression(x[t]),
  xlab = "", ylab = "", typ = "l")
lines(cbind(0:maxT2, c(0, Y1$x[1:maxT2])), lty = 2)
plot(cbind(0:maxT2, c(0, X1[[1]]$ii[1:maxT2])), main = expression(i[t]),
  xlab = "", ylab = "", typ = "l")
lines(cbind(0:maxT2, c(0, Y1$ii[1:maxT2])), lty = 2)
plot(cbind(0:maxT2, c(0, X1[[1]]$ppi[1:maxT2])), main = expression(pi[t]),
  xlab = "", ylab = "", typ = "l", ylim = range(X1[[1]]$ppi[1:maxT2],Y1$ppi[1:maxT2]))
lines(cbind(0:maxT2, c(0, Y1$ppi[1:maxT2])), lty = 2)
plot(cbind(0:maxT2, c(0, -X1[[1]]$ppi[1:maxT2])), main = expression(psi[t]),
  xlab = "", ylab = "", typ = "l", ylim = range(-X1[[1]]$ppi[1:maxT2],Y1$psi[1:maxT2]))
lines(cbind(0:maxT2, c(0, Y1$psi[1:maxT2])), lty = 2)
plot(cbind(0:maxT2, c(0, X1[[1]]$e[1:maxT2])), main = expression(e[t]),
  xlab = "", ylab = "", typ = "l")

dev.copy2eps(file= figpath %+% "nkm_discretion.eps")

