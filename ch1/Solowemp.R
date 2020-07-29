#=======================================
# Chapter 1, Empirical analysis
#   modified on 2019/10/23
#=======================================

library(xlsx)

dset <- read.xlsx("FqReport1.xlsx", 1)

"%+%" <- function(x, y) paste(x, y, sep = "")
# figpath <- "../../figs/"
figpath <- "./"

library(xtable)

#---------------------------------------

alpha <- 1-mean(c((dset$YWH/dset$GDP11NCY)[5:14],(dset$YWH11N/dset$GDP11NCY)[15:39]))

#yset <- seq(7,32,5)
yset <- seq(9,39,5)

lndy <- log(dset$GDP11CY[yset]/dset$GDP11CY[yset-5])/5*100
lndk <- log(dset$KITFA11[yset-1]/dset$KITFA11[yset-1-5])/5*100
lndl <- log(dset$L[yset]/dset$L[yset-5])/5*100

contk <- alpha*lndk
contl <- (1-alpha)*lndl
conta <- lndy-contk-contl

rslt1 <- cbind(lndy, contk, conta, contl)
colnames(rslt1) <- c("GDP¬’·—¦", "Ž‘–{‚ÌŠñ—^“x", "¶ŽY«‚ÌŠñ—^“x","˜J“­‚ÌŠñ—^“x")
rownames(rslt1) <- c("1984-1988","1989-1993","1994-1998","1999-2003","2004-2008","2009-2013","2014-2018")

xtable(rslt1)

#print(rslt1)

#---------------------------------------

s0 <- (dset$I11NCY/dset$GDP11NCY)
s <- mean(s0[30:39])

delta0 <- c()
for(i in 10:39)
  delta0[i] <- (dset$I11CY[i] - (dset$KITFA11[i]-dset$KITFA11[i-1]))/dset$KITFA11[i-1]

delta <- mean(delta0[10:39])

nexo <- log(dset$L[39]/dset$L[29])/10
gexo <- mean(conta[6:7])/(1-alpha)*0.01

c(alpha, delta, gexo, nexo, s)

kstar <- (s/(gexo+nexo+delta))^(1/(1-alpha))

tmax <- 150

Y0 <- dset$GDP11CY[39]
K0 <- dset$KITFA11[38]
L0 <- dset$L[39]
A0 <- (Y0/K0^alpha)^(1/(1-alpha))/L0

A <- c(1, cumprod(rep(1+gexo, tmax-1)))*A0
L <- c(1, cumprod(rep(1+nexo, tmax-1)))*L0

Yest <- K0^alpha*(A[1]*L[1])^(1-alpha)
Kest <- dset$KITFA11[39]

for( i in 2:tmax) {
  Yest[i] <- Kest[i-1]^alpha*(A[i]*L[i])^(1-alpha)
  Kest[i] <- (1-delta)*Kest[i-1] + s*Yest[i]
}

yplot <- 2018:2167
yplot1 <- yplot-2017

Ystar <- kstar^alpha*A*L
kest <- c(K0, Kest[-tmax])/A/L


if (.Platform$OS.type == "windows") windows(10, 5)
par(mfrow = c(1,2))
par(ps = 15)
par(mai = c(0.85, 0.68*0.75, 0.68*0.75, 0.35*1.5)*1)
plot(cbind(yplot,Yest[yplot1]/10^3), typ ="l", xlab = "Year", main = expression(Y[t]), xlim = range(yplot), ylim = range(Ystar/1000), log="y")
lines(cbind(yplot,Ystar[yplot1]/10^3), lty = 2)
axis(side=1, labels = 2018, at = 2018,las = 1)


plot(cbind(yplot,kest[yplot1]), typ ="l", ylim = range(kest), xlab = "Year", main = expression(k[t]), xlim = range(yplot))
lines(cbind(yplot,kstar), lty = 2)
axis(side=1, labels = 2018, at = 2018,las = 1)
axis(side=4, labels = expression("k*"), at = kstar,las = 2)

dev.copy2eps(file= figpath %+% "solforecast1.eps")
dev.off()

