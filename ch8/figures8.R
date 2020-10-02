#=======================================
# Chapter 8, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

#----------------------------
# 主な確率分布

f1 <- function(x) dnorm(x, 0, 1)
f2 <- function(x) dgamma(x, 1, 1)
f3 <- function(x) dbeta(x, 2, 2)
f4 <- function(x) dunif(x, -1, 1)

par(mfrow = c(2,2), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
curve(f1, -4, 4, main = expression(paste("N(",mu," = 0, ", sigma^2," = 1 )")))
abline(h = 0, col = 'grey')

curve(f2, 0, 5, main = expression(paste("Ga(",s," = 1, ", r," = 1 )")),xlim = c(0,4))
abline(h = 0, col = 'grey')
lines(cbind(c(0, 0), c(0, 1)),lty=2)

curve(f3, 0, 1, main = expression(paste("Be(",a," = 2, ", b," = 2 )")))
abline(h = 0, col = 'grey')

curve(f4, -1, 1, main = expression(paste("Unif(",a," = -1, ", b," = 1 )")),ylim = c(0, 1))
lines(cbind(c(-1.0, -1.0), c(0, 0.5)),lty=2)
lines(cbind(c(1.0, 1.0), c(0, 0.5)),lty=2)
abline(h = 0, col = 'grey')

dev.copy2eps(file= figpath %+% "distribution.eps")

#----------------------------
# ベイズ推定の例（ベルヌーイ過程）

N <- c(5, 20, 100, 1000)
theta <- 0.6
set.seed(101)
vec <- runif(N[4])
rslt <- rep(0, N[4])
rslt[vec < theta] <- 1

a0 <- 2
b0 <- 2
x1 <- seq(0,1,length = 200)
y1 <- dbeta(x1, a0, b0)
fb1 <- function(x) dbeta(x, a0+sum(rslt[1:N[1]]), b0+N[1]-sum(rslt[1:N[1]]))
fb2 <- function(x) dbeta(x, a0+sum(rslt[1:N[2]]), b0+N[2]-sum(rslt[1:N[2]]))
fb3 <- function(x) dbeta(x, a0+sum(rslt[1:N[3]]), b0+N[3]-sum(rslt[1:N[3]]))
fb4 <- function(x) dbeta(x, a0+sum(rslt[1:N[4]]), b0+N[4]-sum(rslt[1:N[4]]))

par(mfrow = c(2,2), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

curve(fb1, 0, 1, main = "N = 5", ylim = c(0, 5))
lines(cbind(x1, y1), lty = 2)
lines(cbind(c(theta,theta),c(0, 10^4)),col=3)

curve(fb2, 0, 1, main = "N = 20", ylim = c(0, 5))
lines(cbind(x1, y1), lty = 2)
lines(cbind(c(theta,theta),c(0, 10^4)),col=3)

curve(fb3, 0, 1, main = "N = 100", ylim = c(0, 10))
lines(cbind(x1, y1), lty = 2)
lines(cbind(c(theta,theta),c(0, 10^4)),col=3)

curve(fb4, 0, 1, main = "N = 1000", ylim = c(0, 30))
lines(cbind(x1, y1), lty = 2)
lines(cbind(c(theta,theta),c(0, 10^4)),col=3)

dev.copy2eps(file= figpath %+% "bernoulli.eps")

#----------------------------
# マルコフ過程の例

set.seed(101)
aa <- c(0, -100, 100)
rho <- 0.5

vec1 <- aa[1]
vec2 <- aa[2]
vec3 <- aa[3]

for(i in 1:999) {
  vec1[i+1] <- rho*vec1[i] + rnorm(1)
  vec2[i+1] <- rho*vec2[i] + rnorm(1)
  vec3[i+1] <- rho*vec3[i] + rnorm(1)
}

if (.Platform$OS.type == "windows") windows(10,7)
par(ps = 15)
plot(vec1[1:200], typ = "l", xlab = "", ylab = "", main = expression(a[t]), ylim = c(-5, 5))
lines(vec2[1:200], col = 2,lty=2)
lines(vec3[1:200], col = 4,lty=4)

dev.copy2eps(file= figpath %+% "markov.eps")

#----------------------------
# MHアルゴリズムの例

set.seed(102)

ff <- function(x) log(dgamma(x, 5, 1))
x2 <- seq(0,15,length = 200)
y2 <- dgamma(x2, 5, 1)

niter <- 1200
chain4 <- rgamma(niter, 5, 1)

Metropolis <- function(init, sig2) {
  chain <- init[1]
  for(i in 2:niter) {
    cand <- chain[i-1] + rnorm(1, sd = sqrt(sig2))
    prob <- min(exp(ff(cand)-ff(chain[i-1])), 1)
    if (runif(1) < prob) {
      chain[i] <- cand
    } else {
      chain[i] <- chain[i-1]
    }
  }
  return(chain)
}

chain1 <- Metropolis(9, 0.1)
chain2 <- Metropolis(9, 20)
chain3 <- Metropolis(9, 5000)

rejrate1 <- sum(chain1[201:1200] == chain1[201:1200-1])/1000
rejrate2 <- sum(chain2[201:1200] == chain2[201:1200-1])/1000
rejrate3 <- sum(chain3[201:1200] == chain3[201:1200-1])/1000
1-c(rejrate1,rejrate2,rejrate3)

if (.Platform$OS.type == "windows") windows(7,7)

par(mfrow = c(2,2), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
par(family="Japan1GothicBBB") 

hist(chain1[201:1200], xlab = "", main = expression(paste(sigma^2, " = 0.1")), probability = T, xlim = c(0, 15))
lines(cbind(x2, y2))

hist(chain2[201:1200], xlab = "", main = expression(paste(sigma^2, " = 20")), probability = T, xlim = c(0, 15))
lines(cbind(x2, y2))

hist(chain3[201:1200], xlab = "", main = expression(paste(sigma^2, " = 5000")), probability = T, xlim = c(0, 15))
lines(cbind(x2, y2))

hist(chain4[201:1200], xlab = "", main = "(参考)", probability = T, xlim = c(0, 15), ylim = c(0, 0.2))
lines(cbind(x2, y2))

dev.copy2pdf(file=figpath %+% "metropolis.pdf",family="Japan1GothicBBB")

