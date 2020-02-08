#==================================================
# Chapter 6, RBC model (Chebyshev approximation)
#   modified on 2019/10/19
#==================================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

beta <- 0.99
alpha <- 0.3
gamma <- 1.0
delta <- 0.025
rho <- 0.9
mu <- 1.0

r_ss <- 1/beta-1+delta
k_l <- (r_ss/alpha)^(1/(alpha-1))
w_ss <- (1-alpha)*(k_l)^alpha
c_l <- k_l^alpha-delta*k_l
l_ss <- ((w_ss/mu/(1+gamma)) * c_l^(-1))^(1/(gamma+1))
k_ss <- k_l*l_ss
c_ss <- c_l*l_ss

U <- function (x) log(x[1]) - mu * x[2]^(gamma+1)

g <- function (s, x, epsilon) {
  wt <- (1-alpha) * s[1]^alpha / x[2]^alpha * s[2]
  rt <- alpha * x[2]^(1-alpha) / s[1]^(1-alpha)  * s[2] 
  s1 <- x[2] * wt + s[1] * rt + (1 - delta) * s[1] - x[1]
  s2 <- exp(rho * log(s[2]) + epsilon)
  c(s1, s2)
}

y <- 0
for(i in 0:500)
  y <- y + beta^i * U(c(c_ss,l_ss))

#---------------------------------------
# チェビシェフ近似補間（２次）
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

fn2 <- function(x, cc, nt){
  size <- dim(x)
  lTk <- vector("list", size[2])
  for (i in 1:size[2]){
    Tk <- matrix(1, size[1], 1)
    Tk <- cbind(Tk, x[, i])
    for (k in 2:(nt-1))
     Tk <- cbind(Tk, 2*x[, i]*Tk[, k]-Tk[, k-1])
    lTk[[i]] <- Tk
  }
  Tkk <- c()
  for (i in 1:size[1]) {
    rTk <- lTk[[1]][i,]
    for (k in 2:size[2])
      rTk <- lTk[[k]][i,] %x% rTk
    Tkk <- rbind(Tkk, rTk)
  }
  rownames(Tkk) <- NULL
  Tkk %*% cc
}

xj2 <- function(n = 10){
  as.matrix(expand.grid(xj(n),xj(n)))
}

trans2 <- function(x, a, b) cbind(trans(x[, 1], a[1], b[1]), trans(x[, 2], a[2], b[2]))
rev2 <- function(y, a, b) cbind(rev(y[, 1], a[1], b[1]), rev(y[, 2], a[2], b[2]))

cheb2x <- function(x, coef, a, b, n, vc, z) {
  ans <- fn2(trans2(matrix(x, 1), a, b), coef, n) 
  ans
}

#---------------------------------------
# Policy Iteration
#---------------------------------------

getcoef <- function(nt, f) {
  ngrid <- nt*nt
  at <- c(1.0, 0.90)
  bt <- c(30.0, 1.10)
  
  node <- cbind(rev(xj(nt), at[1], bt[1]), rev(xj(nt), at[2], bt[2]))
  xjgrid <- rev2(xj2(nt), at, bt)
  apxfun <- function(x, coef) fn2(trans2(x, at, bt), coef, nt)
  
  niter <- 1000
  eps <- 0.001
  
  # initialization
  cc <- rep(0, ngrid)
  v <- rep(0, ngrid)
  coef1 <- rep(0, ngrid)
  coef2 <- rep(0, ngrid)
  
  hl <- function(ct, s) 
    ((1-alpha)*s[1]^alpha*s[2]/(gamma+1)/mu/ct)^(1/(alpha+gamma))
  
  for(iter in 1:niter) {
    print(iter)
    v1 <- c(v)
    c1 <- c(cc)
    
    for(i in 1:ngrid) {
      obj <- function(ct) {
        st <- c(xjgrid[i, ])
        lt <- hl(ct, st)
        xt <- c(ct, lt)
        sa <- matrix(g(st,xt , 0), 1)
        - U(xt) - beta*f(sa)
      }
      opt <- optimize(obj, c(0.001, (1 - delta)*xjgrid[i, 1]))
      cc[i] <- opt$minimum
    }
    fobj1 <- function(coef) cc - apxfun(xjgrid, coef)
    coef1 <- nleqslv(coef1, fobj1)$x
    hc <- function(x) cheb2x(x, coef1, at, bt, nt, node, cc)
    
    for(i in 1:ngrid) {
      st <- c(xjgrid[i, ])
      ct <- max(min(hc(st), (1 - delta) * st[1]), 0.1)
      xt <- c(ct, hl(ct, st))
      y <- U(xt)
      for(k in 1:500) {
        st <- g(st, xt, 0)
        ct <- max(min(hc(st), (1 - delta) * st[1]), 0.1)
        xt <- c(ct, hl(ct, st))
        y <- y + beta^k * U(xt)
      }
      v[i] <- y
    }
    if (identical(abs(c(cc) - c1) < eps ,rep(T, nt*nt))) break
    
    fobj2 <- function(coef) v - apxfun(xjgrid, coef)
    coef2 <- nleqslv(coef2, fobj2)$x
    f <- function(x) cheb2x(x, coef2, at, bt, nt, node, v)
    
    persp(node[, 1], node[, 2], matrix(v, nt, nt))
    print(abs(c(cc) - c1))
  }
  return(list(coef1, coef2, cc))
}

f <- function (x) -45 + sqrt(x[1]) + sqrt(x[2])
at <- c(1.0, 0.90)
bt <- c(30.0, 1.10)
nt1 <- 10
coefset1 <- getcoef(10, f)

f2 <- function(x) cheb2x(x, coefset1[[2]], at, bt, nt1)
nt2 <- 15

coefset2 <- getcoef(nt2, f2)
hc <- function(x) cheb2x(x, coefset2[[1]], at, bt, nt2)
hl <- function(ct, s) 
    ((1-alpha)*s[1]^alpha*s[2]/(gamma+1)/mu/ct)^(1/(alpha+gamma))
cc <- coefset2[[3]]

#---------------------------------------
# Plot
#---------------------------------------

nt <- nt2
node <- cbind(rev(xj(nt), at[1], bt[1]), rev(xj(nt), at[2], bt[2]))

st <- c(k_ss, 1.0)
kt <- c(k_ss)
for (t in 1:10000) {
  ct <- hc(st)
  lt <- hl(ct, st)
  st <- g(st, c(ct, lt), 0)
  kt[t+1] <- st[1]
}
k_sn <- kt[t+1]

nsim <- 101
X <- matrix(NA, nsim, 4)
X[1, 3:4] <- c(k_sn, 1.0)
X[1, 1] <- hc(X[1, 3:4])
X[1, 2] <- hl(X[1, 1], X[1, 3:4])
X[2, 3:4] <- g(X[1, 3:4], X[1, 1:2], 0.01)

for (i in 2:nsim) {
  X[i, 1] <- hc(X[i, 3:4])
  X[i, 2] <- hl(X[i, 1], X[i, 3:4])
  if (i < nsim)
    X[i+1, 3:4] <- g(X[i, 3:4], X[i, 1:2], 0)
}

windows(7, 7)
par(ps = 12, mai = c(0.85, 0.68, 0.68, 0.35)*0.3)
persp(node[, 1], node[, 2], matrix(cc, nt, nt),
      theta = -30, phi = 30, expand = 0.5, ticktype = "detailed",
      col = "lightgray", border = "red", xlab = "K", ylab = "A", zlab = "C")

dev.copy2eps(file= figpath %+% "polit_rbc_pol.eps")

par(mfrow=c(2,2), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
plot(cbind(0:(nsim-1), (log(X[, 1])-log(c_ss))*100), typ = "l", col = 2, lwd = 1, xlab = "", ylab = "", main = expression(C[t]))
plot(cbind(0:(nsim-1), (log(X[, 3])-log(k_ss))*100), typ = "l", col = 2, lwd = 1, xlab = "", ylab = "", main = expression(K[t]))
plot(cbind(0:(nsim-1), (log(X[, 2])-log(l_ss))*100), typ = "l", col = 2, lwd = 1, xlab = "", ylab = "", main = expression(L[t]))
plot(cbind(0:(nsim-1), log(X[, 4])*100), typ = "l", col = 2, lwd = 1, xlab = "", ylab = "", main = expression(A[t]))
dev.copy2eps(file= figpath %+% "polit_rbc_ckla.eps")

