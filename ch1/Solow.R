#=======================================
# Chapter 1, Solow model
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

#--------------------
# ƒ\ƒ[ƒ‚ƒfƒ‹
#--------------------

alpha <- 0.30
delta <- 0.1
s <- 0.2
g <- 0.02
n <- -0.01

fz1 <- function(k) s*k^alpha
fz2 <- function(k) (g + n + delta)*k
dk <- function(k) (s*k^alpha-(g + n + delta)*k)/(1+g)/(1+n)

ks <- (s/(g+n+delta))^(1/(1-alpha))

kk <- seq(0,10,length = 201)

par(ps = 15)
par(mai = c(0.85, 0.85, 0.35, 0.35))
plot(cbind(kk, fz1(kk)), typ = "l", xlim = c(0, 5.7), ylim = c(0,0.6),
     xlab = expression(k[t]), ylab = expression(z))
lines(cbind(kk, fz2(kk)))
abline(v = ks, lty = 2)
abline(h = 0)

text(1, 0.26, labels = expression(paste("z = s ",k[t]^alpha)))
text(4.8, 0.42  , labels = expression(paste("z = ",(g+n+delta)," ", k[t])))
text(2.9, 0.5, labels = expression(paste(k[t]," = ","k*")))

axis(side=1, labels = expression(k[1]), at = ks*.5)

k1 <- ks*.5
segments(k1, 0, ks*.5, fz2(k1), lty = 2)
arrows(k1, fz2(k1),k1, fz1(k1),length = 0.1,col=4,lwd = 2)
arrows(k1, fz1(k1),k1+(fz1(k1)-fz2(k1))*6 ,fz1(k1),length = 0.1,col=4,lty=2,lwd = 2)

dev.copy2eps(file= figpath %+% "solow1.eps")

#-----------------------------

ka  <- ks*.5
dka <- NA
for(i in 2:110) {
  dka[i] <- dk(ka[i-1])
  ka[i] <- ka[i-1] + dka[i]
}

par(ps = 17.5)
par(mfrow = c(2, 1))
par(mai = c(0.85*0.5, 0.68*0.75, 0.68*0.75, 0.5)*1.1)
plot(ka, xlim = c(0, 100), xlab = "", main = expression(k[t]), ylab = "", typ = "l")
abline(h = k1, lty = 2, col = 2)
axis(side=4, labels = expression(k[1]), at = k1, las = 2)
axis(side=4, labels = expression("k*"), at = ks, las = 2)

plot(dka, xlim = c(0, 100), xlab = "", main = expression(paste(Delta, k[t])), ylab = "", typ = "l")

dev.copy2eps(file= figpath %+% "solow1a.eps")

#-----------------------------

k2 <- ks*2
kb <- k2
dkb <- NA
for(i in 2:110) {
  dkb[i] <- dk(kb[i-1])
  kb[i] <- kb[i-1] + dkb[i]
}

plot(kb, xlim = c(0, 100), xlab = "", main = expression(k[t]), ylab = "", typ = "l")
abline(h = k2, lty = 2, col = 2)
axis(side=4, labels = expression(k[1]), at = k2, las = 2)
axis(side=4, labels = expression("k*"), at = ks, las = 2)

plot(dkb, xlim = c(0, 100), xlab = "", main = expression(paste(Delta, k[t])), ylab = "", typ = "l")

dev.copy2eps(file= figpath %+% "solow1b.eps")

#-----------------------------


fs <- function(s) (s^(alpha/(1-alpha))-s^(1/(1-alpha)))/(g+n+delta)^(alpha/(1-alpha))

par(ps = 15)
par(mai = c(0.85, 0.85, 0.35, 0.35))
par(mfrow = c(1, 1))
curve(fs,0,1,ylim=c(0, 1.5), xlab = "s", ylab = expression("c*"))
abline(h = 0, v = 0)

cg <- fs(alpha)
lines(cbind(c(0,1),c(cg,cg)),lty = 2)

points(alpha,cg, pch=19)
text(alpha+0.075,cg+0.075, labels = expression(paste("(",s[g],", ",c[g],")")))

text(0.8, 0.65, labels = expression("c* = f(s)"))
text(0.85, 1.0, labels = expression(paste("c* = f(",s[g],")")))

dev.copy2eps(file= figpath %+% "solow2.eps")

#-----------------------------

fct <- function(k, s) (1-s)*k^alpha/(1+g)/(1+n)
fks <- function(s) (s/(g+n+delta))^(1/(1-alpha))
dks <- function(k, s) (s*k^alpha-(g + n + delta)*k)/(1+g)/(1+n)

sA0 <- 0.2
sA1 <- alpha

ktA <- fks(sA0)
ctA <- fct(ktA[1], sA0)

for(i in 2:120)
  ktA[i] <- ktA[i-1] + dks(ktA[i-1], sA1)

ctA[2:120] <- fct(ktA[-1], sA1)


par(ps = 17.5)
par(mfrow = c(2, 1))
par(mai = c(0.85*0.5, 0.68*0.75, 0.68*0.75, 0.5)*1.1)
plot(ctA, xlim = c(0, 100), xlab = "", main = expression(c[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(c[0]), at = ctA[1], las = 2)

plot(ktA, xlim = c(0, 100), xlab = "", main = expression(k[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(k[0]), at = ktA[1], las = 2)

dev.copy2eps(file= figpath %+% "solcaseA.eps")

#-----------------------------

sB0 <- alpha
sB1 <- 0.4

ktB <- fks(sB0)
ctB <- fct(ktB[1], sB0)

for(i in 2:120)
  ktB[i] <- ktB[i-1] + dks(ktB[i-1], sB1)

ctB[2:120] <- fct(ktB[-1], sB1)

plot(ctB, xlim = c(0, 100), xlab = "", main = expression(c[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(c[0]), at = ctB[1], las = 2)

plot(ktB, xlim = c(0, 100), xlab = "", main = expression(k[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(k[0]), at = ktB[1], las = 2)

dev.copy2eps(file= figpath %+% "solcaseB.eps")

#-----------------------------

sC0 <- 0.4
sC1 <- alpha

ktC <- fks(sC0)
ctC <- fct(ktC[1], sC0)

for(i in 2:120)
  ktC[i] <- ktC[i-1] + dks(ktC[i-1], sC1)

ctC[2:120] <- fct(ktC[-1], sC1)

plot(ctC, xlim = c(0, 100), xlab = "", main = expression(c[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(c[0]), at = ctC[1], las = 2)

plot(ktC, xlim = c(0, 100), xlab = "", main = expression(k[t]), ylab = "", typ = "l")
axis(side=4, labels = expression(k[0]), at = ktC[1], las = 2)


dev.copy2eps(file= figpath %+% "solcaseC.eps")

