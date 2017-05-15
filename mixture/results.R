#Compare results

library(INLA)
library(rjags)


#Load results
load("oldfaith.RData")
load("jags-gaussian-ngrp-2.RData")


#Display mean of the labels
plot(apply(res.jags$Z, 1, mean), apply(zz, 2, mean), 
  xlab = "MCMC", ylab = "INLA within MCMC"
)
abline(0, 1)


#Marginals of the fixed effects and hyperparameters
margs.fixed <- lapply(inlamh.res$model.sim, function(X) {
  X$model$marginals.fixed})

margs.hyper <- lapply(inlamh.res$model.sim, function(X) {
  X$model$marginals.hyperpar})

n.margs <- length(margs.fixed)

ws <- rep(1/n.margs, n.margs)

get.margs <- function(margs, ws) {
  lapply(1:length(margs[[1]]), function(X) {
   aux <- lapply(margs, function(Y){ Y[[X]]})
  INLABMA:::fitmargBMA(aux, ws)
})
}

margs.fixed <- get.margs(margs.fixed, ws)
margs.hyper <- get.margs(margs.hyper, ws)


#Probabilites
probs.mcmc <- apply(res.jags$Z == 1, 1, mean)
probs.inlamcmc <- apply(zz == 1, 2, mean)

idx <- order (yy)

pdf(file = "class.pdf")

plot(yy[idx], probs.inlamcmc[idx], main = "Probability of being in group 1",
  xlab = "Eruption time", ylab = "Probability", type = "l")

lines(yy[idx], probs.mcmc [idx], lty = 2)
legend("topright", lty = 1:2, legend = c("INLA wihin MCMC", "MCMC"), bty = "n",
   cex = 1)

dev.off()


pdf(file = "mixtures.pdf")

par(mfrow = c(2, 2))

for(i in 1:2) {
plot(margs.fixed[[i]], type = "l", xlab = "", ylab = "",
  main = substitute(mu[i], list (i = i)))
lines(density(res.jags$mu[i,,]), lty = 2)
legend("topright", lty = 1:2, legend = c("INLA w/ MCMC", "MCMC"), bty = "n",
   cex = 0.8)

plot(margs.hyper[[i]], type = "l",  xlab = "", ylab = "",
  main = substitute(tau[i], list (i = i)))
lines(density(res.jags$prec[i,,]), lty = 2)
legend("topright", lty = 1:2, legend = c("INLA within MCMC", "MCMC"), bty = "n",
  cex = 0.8)

}

dev.off()

