#=======================================
# Chapter 4, Drawing data
#   modified on 2019/10/19
#=======================================

"%+%" <- function(x, y) paste(x, y, sep = "")
figpath <- "../../figs/"

dat <- read.csv("FqReport4.csv")

# peak   bottom
# 1980Q1 1983Q1
# 1985Q2 1986Q4
# 1991Q1 1993Q4
# 1997Q2 1999Q1
# 2000Q4 2002Q3
# 2008Q1 2009Q1
# 2012Q1 2002Q4

peak <- c(1980, 1985.25, 1991, 1997.25, 2000.75, 2008, 2012)
bottom <- c(1983, 1986.75, 1993.75, 1999, 2002, 2009, 2012.75)

cpeak0 <- c("1980Q1", "1985Q2", "1991Q1", "1997Q2", "2000Q4", "2008Q1", "2012Q1")
cbottom0 <- c("1983Q1", "1986Q4", "1993Q4", "1999Q1", "2002Q1", "2009Q1", "2012Q4")

cpeak <- c("80Q1", "85Q2", "91Q1", "97Q2", "2000Q4", "08Q1", "12Q1")
cbottom <- c("83Q1", "86Q4", "93Q4", "99Q1", "02Q1", "09Q1", "12Q4")

gap.ts <- ts(dat$GDPGAP, start = 1980, frequency = 4)
roh.ts <- ts(dat$ROH, start = 1980, frequency = 4)

#---------

windows(7, 4)
par(ps = 13)
par(mai = c(0.6, 0.85, 0.6, 0.35))
plot(gap.ts, main = "", ylab = "%", xlab = "")
abline(v = bottom, lty = 2, col = grey(0.5))
abline(v = peak, col = grey(0.5))

par(ps = 13)
axis(side=3, labels = cpeak, at = peak)

dev.copy2eps(file= figpath %+% "gap.eps")

#---------

par(ps = 13)
plot(roh.ts, main = "", ylab = "2015 = 100", xlab = "")
abline(v = bottom, lty = 2, col = grey(0.5))
abline(v = peak, col = grey(0.5))

par(ps = 13)
axis(side=3, labels = cpeak, at = peak)

dev.copy2eps(file= figpath %+% "roh.eps")

