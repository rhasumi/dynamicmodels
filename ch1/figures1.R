#=======================================
# Chapter 1, Other figures
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

#--------------------
# 指数関数、対数関数
#--------------------

par(ps = 20)
par(mai = c(0.85, 0.88, 0.35, 0.35))
curve(exp, -3, 3)
dev.copy2eps(file= figpath %+% "exp.eps")

curve(log, -0, 6, ylab = "ln(x)")
dev.copy2eps(file= figpath %+% "log.eps")

#--------------------
# 自然対数の意味
#--------------------
yt <- cumprod(c(1, rep(1.05, 100)))*100
plot(cbind(0:40, yt[1:41]), typ = "l", xlab = "t", ylab = expression(Y[t]))
dev.copy2eps(file= figpath %+% "exln1.eps")

plot(cbind(0:40, log(yt[1:41])), typ = "l", xlab = "t", ylab = expression(paste(paste("ln(", Y[t]), ")")))
dev.copy2eps(file= figpath %+% "exln2.eps")

#--------------------
# 生産関数
#--------------------
par(ps = 20)
fk <- function(K) K^0.3
curve(fk, 0, 3, xlab = "K", ylab = "f(K)")
dev.copy2eps(file= figpath %+% "fk1.eps")
dev.copy2pdf(file= figpath %+% "fk1.pdf")

curve(fk, 0, 3, xlab = "K", ylab = "f(K)")
points(cbind(c(1,1.5,2),fk(c(1,1.5,2))), pch=19)

#--------------------
# 微分の意味
#--------------------
b1 <- fk(2)-1
b2 <- (fk(1.5)-1)/0.5
# b3 <- (fk(1.1)-1)/0.1
b4 <- 0.3

abline(1-b1, b1, lty = 2, col = 4)
abline(1-b2, b2, lty = 3, col = 4)
abline(1-b4, b4, col = 4)

legend(1.7, 0.3, legend = c("傾き = 0.23    ","傾き = 0.26"    ,"傾き = 0.3    "), lty = c(2,3,1), col = c(4,4,4))
dev.copy2eps(file= figpath %+% "fk2.eps", family="Japan1GothicBBB")
dev.copy2pdf(file= figpath %+% "fk2.pdf", family="Japan1GothicBBB")


par(ps = 15)
fr <- function(r) log(1+r)
curve(fr, -1, 2, xlab = "r", ylab = "y")
abline(0, 1, lty = 1, col = 4)

text(-0.5, -2.5, "y = ln(1+r)")
text(-0.75, -0.5, "y = r")

dev.copy2eps(file= figpath %+% "ln1r.eps")

