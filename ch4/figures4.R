#=======================================
# Chapter 4, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

A0 <- 1.0
A1 <- 1.01
rhos <- c(0, 0.5, 0.9, 1.0)
Len <- 100

ar1 <- function(rho, A1, len) {
  a <- log(A1)
  for(i in 2:len)
    a[i] <- rho*a[i-1]
  exp(a)
}

ans1 <- c(A0, ar1(rhos[1], A1, Len))
ans2 <- c(A0, ar1(rhos[2], A1, Len))
ans3 <- c(A0, ar1(rhos[3], A1, Len))
ans4 <- c(A0, ar1(rhos[4], A1, Len))

windows(7, 7)
par(mfrow = c(2,2))
par(ps = 15)
par(mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

plot(cbind(0:Len, ans1), xlim = c(-5, 50), typ = "l", xlab = "t", ylab = expression(A[t]), main = expression(paste(rho, " = 0")))
plot(cbind(0:Len, ans2), xlim = c(-5, 50), typ = "l", xlab = "t", ylab = expression(A[t]), main = expression(paste(rho, " = 0.5")))
plot(cbind(0:Len, ans3), xlim = c(-5, 50), typ = "l", xlab = "t", ylab = expression(A[t]), main = expression(paste(rho, " = 0.9")))
plot(cbind(0:Len, ans4), xlim = c(-5, 50), typ = "l", xlab = "t", ylab = expression(A[t]), main = expression(paste(rho, " = 1.0")))

dev.copy2eps(file= figpath %+% "rbcar.eps")

