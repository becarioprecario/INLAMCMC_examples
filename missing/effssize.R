#Test effective sample size

library(INLA)
library(rjags)
library(coda)

#Load data
load("INLA-MH-missing.RData")
load("jags.RData")


#Create mcmc object from INLA+MCMC
b.sim.mcmc <- as.mcmc(x.sim)

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
jags.mcmc <- as.mcmc(t(jm1.samp$x[idx.mis,,]))

effsize.jags <- get.ess(jags.mcmc)
 

#Display results
pdf(file = "effssize-missing.pdf")

plot(n.iter, effsize.inlamcmc, type = "l", 
  xlab = "Number of iterations",
  ylab = "Effective sample size",
  main = "Minimum effective sample size (linear regression)")
lines(n.iter, effsize.jags, lty = 2)

legend("topleft", legend = c("INLA within MCMC", "MCMC"),
  lty = 1:2, bty = "n")

dev.off()

