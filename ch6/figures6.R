#=======================================
# Chapter 6, Dynamic programming
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

xj <- function(n = 9){
  -cos(0.5*(2*(1:n)-1)*pi/n)
}

phi <- function(x, j) {
  if (j == 0) {
    y <- rep(1, length(x))
  } else if (j == 1) {
    y <- x
  } else {
    y <- 2*x*phi(x, j-1) - phi(x, j-2)
  }
  return(y)
}
nodes <- seq(-1,1,length = 100)

windows(7,2.5)
par(ps = 10, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)
plot(cbind(xj(),rep(0,9)),main = expression(z["i"]), xlab = "", ylab = "",axes = F)
axis(1)
axis(2, labels = "", at = 0)
box()
dev.copy2eps(file= figpath %+% "chebnode.eps")

windows(7,7)
par(mfrow = c(3,3), ps = 15, mai = c(0.85, 0.68, 0.68, 0.35)*0.5)

ff <- function(x) phi(x, j=0)
curve(ff, -1, 1, main = expression(paste(phi[0],"(z)")))

ff <- function(x) phi(x, j=1)
curve(ff, -1, 1, main = expression(paste(phi[1],"(z)")))

ff <- function(x) phi(x, j=2)
curve(ff, -1, 1, main = expression(paste(phi[2],"(z)")))

ff <- function(x) phi(x, j=3)
curve(ff, -1, 1, main = expression(paste(phi[3],"(z)")))

ff <- function(x) phi(x, j=4)
curve(ff, -1, 1, main = expression(paste(phi[4],"(z)")))

ff <- function(x) phi(x, j=5)
curve(ff, -1, 1, main = expression(paste(phi[5],"(z)")))

ff <- function(x) phi(x, j=6)
curve(ff, -1, 1, main = expression(paste(phi[6],"(z)")))

ff <- function(x) phi(x, j=7)
curve(ff, -1, 1, main = expression(paste(phi[7],"(z)")))

ff <- function(x) phi(x, j=8)
curve(ff, -1, 1, main = expression(paste(phi[8],"(z)")))

dev.copy2eps(file= figpath %+% "chebfun.eps")

