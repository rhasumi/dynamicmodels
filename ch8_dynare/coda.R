#=======================================
# Chapter 8, Post-processing
#   modified on 2019/10/19
#=======================================

library(coda)
library(xtable)

chain1 <- read.csv("chain1.csv")
chain2 <- read.csv("chain2.csv")

colnames(chain1) <- c("epsilon", "z", "e", "gamma", "varrho", "phi_pi", "phi_y", "rho_A", "rho_v")
colnames(chain2) <- c("epsilon", "z", "e", "gamma", "varrho", "phi_pi", "phi_y", "rho_A", "rho_v")

vorder <- c(7, 8, 9, 1, 2, 3, 4, 5)
varnames <- c("$gamma$", "$varrho$", "$phi_pi$", "$phi_y$", "$rho_A$", "$rho_v$", "$sigma_varepsilon$", "$sigma_z$","$sigma_e$")

rslt <- mcmc.list(mcmc(chain1[25001:125000, ]), mcmc(chain2[25001:125000, ]))
sum <- summary(rslt)
gew <- geweke.diag(rslt)

ans0 <- as.matrix(sum$statistics[,1:2])
ans0 <- cbind(ans0, gew[[1]]$z, gew[[2]]$z)
ans <- ans0[c(4:9, 1:3), ]

rownames(ans) <- varnames
colnames(ans) <- c("Mean","SD","Z(1)", "Z(2)")

print(ans)

cbind(1-rejectionRate(rslt[[1]][,1]), 1-rejectionRate(rslt[[2]][,1]))

xtable(ans,digits = 3)
