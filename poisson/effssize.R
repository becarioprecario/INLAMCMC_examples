#Test effective sample size

library(INLA)
library(rjags)
library(coda)

#Load data
load("INLA-MHbivpo.RData")
load("jags.RData")


#Create mcmc object from INLA+MCMC
b.sim.mcmc <- as.mcmc(b.sim)

#Compute eff. s.size for different number of iterations
n.iter <- seq(100, 10000, by = 100)

#INLA+MCMC
#samples: mcmc object

get.ess <- function(samples) {
  unlist(lapply(n.iter, function(X) {
    min(effectiveSize(samples[1:X, ]))
  }))
}

effsize.inlamcmc <- get.ess(b.sim.mcmc)

#JAGS
jags.mcmc <- do.call(cbind, 
  lapply(c("alpha", "beta1", "beta2"), function(X) {
      as.vector(jm2.samp[[X]])
  })
)
effsize.jags <- get.ess(jags.mcmc)
effsize.jags2 <- get.ess(jags.mcmc[, 2:3])

 

#Display results
pdf(file = "effssize-poisson.pdf")

plot(n.iter, effsize.inlamcmc, type = "l", xlab = "Number of iterations",
  ylab = "Effective sample size",
  main = "Minimum effective sample size (Poisson regression)")
lines(n.iter, effsize.jags, lty = 2)
lines(n.iter, effsize.jags2, lty = 3)


legend("topleft", legend = c(
  expression(paste("INLA within MCMC (",  beta[1], ", ", beta[2], ")")),
 expression(paste("MCMC (", beta[1], ", ", beta[2], ", ", alpha, ")")),
  expression(paste("MCMC (", beta[1], ", ", beta[2], ")"))), lty = 1:3,
  bty = "n")


dev.off()

