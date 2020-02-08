#==================================================
# Chapter 6, RBC model (linear approximation)
#   modified on 2019/10/19
#==================================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

alpha <- 0.3
beta <- 0.99
delta <- 0.025
mu <- 1.0
gamma <- 1.0
rho <- 0.9

Astar <- 1
rstar <- 1/beta + delta - 1
K_L <- (rstar/alpha/Astar)^(1/(alpha-1))
Y_L <- Astar*K_L^alpha
C_L <- Y_L-delta*K_L
wstar <- (1-alpha)*Astar*K_L^alpha
Lstar <- (wstar/(gamma+1)/mu)^(1/(gamma+1))*C_L^(-1/(gamma+1))
Kstar <- K_L*Lstar
Ystar <- Y_L*Lstar
Cstar <- C_L*Lstar
Rstar <- rstar+1

nvar <- 7
Bmat <- matrix(0, nvar, nvar)
Cmat <- matrix(0, nvar, nvar)

Bmat[2,1] <- 1.0
Bmat[2,5] <- -Rstar*beta
Bmat[6,6] <- Kstar
Bmat[7,7] <- 1.0

Cmat[1,1] <- 1.0
Cmat[1,2] <- gamma
Cmat[1,4] <- -1.0
Cmat[2,1] <- 1.0
Cmat[3,2] <- -alpha+1.0
Cmat[3,3] <- -1.0
Cmat[3,6] <- alpha
Cmat[3,7] <- 1.0
Cmat[4,2] <- -alpha
Cmat[4,4] <- -1.0
Cmat[4,6] <- alpha
Cmat[4,7] <- 1.0
Cmat[5,2] <- -alpha+1.0
Cmat[5,5] <- -Rstar/(Rstar-1.0)
Cmat[5,6] <- alpha-1.0
Cmat[5,7] <- 1.0
Cmat[6,1] <- -Cstar
Cmat[6,3] <- Ystar
Cmat[6,6] <- -Kstar*(delta-1.0)
Cmat[7,7] <- rho

A <- solve(Cmat) %*% Bmat

# Policy/Transition Function ‚ð‹‚ß‚éŠÖ”
# ‰Á“¡—Á[2007]‚àŽQÆ

W <- eigen(A)[[2]]
theta <- eigen(A)[[1]]
Q <- solve(W)

ntheta <- length(theta)

## Extract unstable vectors
sjjw <- which(theta - 1.0 < 0)
n <- length(sjjw)
UQ <- Q[sjjw,]

## Extract stable vectors
sjw <- setdiff(c(1:ntheta),sjjw) 
k <- length(sjw)
SQ <- Q[sjw,]

## Extract stable roots
VLL <- theta[sjw]
VL <- diag(1.0/VLL)

# Elements in Q
PA <- UQ[1:n,1:n]
PB <- UQ[1:n,(n+1):(n+k)]
PC <- SQ[1:k,1:n]
PD <- SQ[1:k,(n+1):(n+k)]
P <- -solve(PA) %*% PB
PE <- PC %*% P + PD 

# Solution
AA  <- solve(PE) %*% VL %*% PE

# funcition for simulation
rbc_sim <- function(AA, P, len_t, S_0){
  Ss <- S_0
  S <- matrix(0, nrow=len_t-1, ncol=ncol(AA))
  for (i in 1:(len_t-1)){
    q <- AA %*% Ss
    S[i,] <- t(q)
    Ss <- S[i,]
  }
  SY <- rbind(t(S1), S)
  X <- t(P %*% t(SY))
  return(list(X, SY))
}

# Time span
nsim <- 100

# Initial Values
S1 <- matrix(c(0, 1), ncol=1)

# SIMULATION
simrslt <- rbc_sim(AA, P, nsim, S1) 
X <- simrslt[[1]]
SY  <- simrslt[[2]]

# C,L,Y,w,R,K,A
par(mfrow=c(2,2), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
plot(cbind(0:nsim, c(0, X[, 1])),  typ = "l", col = 4, lwd = 1, xlab = "", ylab = "", main = expression(C[t]))
plot(cbind(0:nsim, c(0, SY[, 1])), typ = "l", col = 4, lwd = 1, xlab = "", ylab = "", main = expression(K[t]))
plot(cbind(0:nsim, c(0, X[, 2])),  typ = "l", col = 4, lwd = 1, xlab = "", ylab = "", main = expression(L[t]))
plot(cbind(0:nsim, c(0, SY[, 2])), typ = "l", col = 4, lwd = 1, xlab = "", ylab = "", main = expression(A[t]))
dev.copy2eps(file= figpath %+% "bk_rbc_ckla.eps")

