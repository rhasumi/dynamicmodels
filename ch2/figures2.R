#=======================================
# Chapter 2, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

#--------------------
# Œø—pŠÖ”
#--------------------
U <- function(a, b) sqrt(a*b)
Ub <- function(a, U) U^2/a

feasible <- c()
amax <- 6
bmax <- 12
budgit <- 600
pa <- 100
pb <- 50
for (a in 0:amax)
  for(b in 0:bmax)
    if (a*pa + b*pb <= budgit)
      feasible <- rbind(feasible, c(a,b))

par(ps = 15)
plot(feasible, xlab = "a", ylab = "b", xlim = c(0,bmax), ylim = c(0,bmax), pch = 20)

dev.copy2eps(file= figpath %+% "feasible1.eps")

Ub1 <- function(a) Ub(a, U(4, 4))
curve(Ub1,0,10, xlab = "a", ylab = "b", xlim = c(0,bmax), ylim = c(0,bmax))
points(feasible, pch = 20)

text(4.3, 4.3, labels = "C")
text(2.3, 8.3, labels = "D")
text(3.3, 6.3, labels = "E")

dev.copy2eps(file= figpath %+% "feasible2.eps")


par(ps = 20)
par(mai = c(0.85, 0.88, 0.35, 0.35))
funsqx <- function(x) sqrt(x)
par(ps = 20, mai = c(0.85, 0.9, 0.68, 0.35))
curve(funsqx, 0, 6, main = "", ylab = expression(sqrt(x)))
dev.copy2eps(file= figpath %+% "sqrtx.eps")

funx2 <- function(x) x^2
par(ps = 20, mai = c(0.85, 1.0, 0.68, 0.35))
curve(funx2, -3, 3, main = "", ylab = expression(x^2))
dev.copy2eps(file= figpath %+% "x2.eps")



#--------------------
# Œø—pÅ‘å‰»
#--------------------

U1 <- function(a, b, U) exp(U-log(a))
U2 <- function(a, b) U1(a, b, log(125/4))

par(ps = 15)
par(mai = c(0.85, 0.85, 0.2, 0.2))
curve(U2, 0, 15, xlim=c(0,13),ylim=c(0,13),xlab = "a", ylab = "b")
lines(cbind(c(0, 0),c(0,100)))
lines(cbind(c(0, 100),c(0,0)))
lines(cbind(c(10, 0),c(0,12.5)))

points(5,6.25, pch=19)
text(7.5,6.25, "(a = 5, b = 6.25)")

text(5, 12, "ln(a) + ln(b) = U*")
#text(10, 4, "U* = ln(a) + ln(b)")
text(11, 1.9, "100a + 80b = 1000")

dev.copy2eps(file= figpath %+% "util1.eps")

