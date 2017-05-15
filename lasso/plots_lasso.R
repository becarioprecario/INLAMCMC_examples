#
#Load data and plot 
#


library(INLA)
library(rjags)
library(glmnet)


#Load data
load("lasso.RData")
load("INLA-lasso.RData")


#Summary stats for coefficients
apply(b.sim, 2, mean)
apply(b.sim, 2, sd)


#Plot with different estimates
pdf(file = "lasso.pdf")

#Bandwidths
bws <- c(0.01, 0.05, 0.01, 0.02, 0.03)

pos <- rep("topleft", 5)
pos[2] <- "topright"
pos[4] <- "topright"
pos[5] <- "topright"

par(mfrow = c(3, 2))
for(i in 1:5) {
  #INLA-MCMC estimates
  plot(density(b.sim[, i], bw = bws[i]), main = colnames(x)[i], xlab = "")
  #MCMC estimates
  lines(density(smp.jags$b[i, , ], bw = bws[i]), lty = 2)
  #Lasso estimates
  abline(v = (as.vector(lasso.coef)[-1])[i], lty = 3)

legend(pos[i], lty = 1:3, bty = "n",
  legend = c("INLA w/ MCMC", "MCMC", "Lasso"))
}

##INLA-MCMC
#plot(density(1/smp.jags$tau[, , 1]), type = "n", main = "Residual variance")
##MCMC estimates
#lines(density(1/smp.jags$tau[, , 1]))
##Lasso estimates (??)
#abline(v = var(y - lasso.fitted), lty = 2)
#legend("topleft", lty = c(1, 1, 2), bty = "n",
#  legend = c("INLA", "MCMC", "Lasso"), col = c("red", "black", "black"))
#
dev.off()


