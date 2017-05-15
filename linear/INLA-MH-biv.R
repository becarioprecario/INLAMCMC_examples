#Implementation of INLA within M-H for Bayesian Inference
#Bivariate inference

#Load libraries
library(INLA)
#INLA:::inla.dynload.workaround()

#INLAMH()
library(INLABMA)


#Simulate data
n <- 100
set.seed(123)
x1 <- runif(n)
x2 <- runif(n)
err <- rnorm(n)
y <- 3 + 2*x1-2*x2+err

d <- data.frame(y = y, x1 = x1, x2 = x2)

save(file="databiv.RData", list = c("d"))

#INLA models
m1 <- inla(y ~ 1+x1+x2, data = d)

#Fit linear model with R-INLA with a fixed beta
#b Vector of length 2
fit.inla <- function(data, b) {

   data$oset <- b[1] * data$x1 +b[2]*data$x2

   res <- inla(y ~1+offset(oset), data = data) 

   return(list(mlik = res$mlik[1,1], model = res))
}

#Test model fitting
fit.inla(d, c(1,1))


#Proposal x -> y
#density
dq.beta <- function(x, y, sigma = .75, log =TRUE) {
	sum(dnorm(y, mean = x, sd = sigma, log = log))
}
#random
rq.beta <- function(x, sigma = .75) {
	rnorm(length(x), mean = x, sd = sigma)
}

#Prior for beta
prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
   sum(dnorm(x, mean = 0, sd= sigma, log = log))
}


#Run simulations
inlamh.res <- INLAMH(d, fit.inla, c(0, 0), rq.beta, dq.beta, prior.beta,
  n.sim = 10000, n.burnin = 500, n.thin = 10)


#Show results
b.sim <- do.call(rbind, inlamh.res$b.sim)
model.sim <- inlamh.res$model.sim

save(file = "INLA-MHbiv.RData", list = ls())

#Load JAGS results
#
#Run source("jags.R", echo = TRUE) first...
#
library(rjags)
load("jags.RData")

#Compute conturs
library(MASS)
z.inla <- kde2d(b.sim[, 1], b.sim[, 2])
z.mcmc <- kde2d(jm2.samp$beta1[1, , ], jm2.samp$beta2[1, , ])


pdf(file = "INLA-MHbiv.pdf", width = 10, height = 5)
par(mfrow = c(1, 3))
plot(b.sim)
#abline(h = 2, col = "red")

#Density beta1
plot(density(b.sim[, 1]))
lines(m1$marginals.fixed[[2]], col ="red")
legend("topleft", lty = 1, legend = c("INLA", "MCMC"), col = c("red", "black"))

#Density beta2
plot(density(b.sim[, 2]))
lines(m1$marginals.fixed[[3]], col ="red")
legend("topleft", lty = 1, legend = c("INLA", "MCMC"), col = c("red", "black"))

dev.off()


print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 

#Summary statistics
print("Summary statistics of beta:")
c(mean(b.sim[, 1]), sd(b.sim[, 1]), quantile(b.sim[, 1], c(0.025, .5, .975)) )
c(mean(b.sim[, 2]), sd(b.sim[, 2]), quantile(b.sim[, 2], c(0.025, .5, .975)) )
m1$summary.fixed[-1,]

#We need to add the prior of the paramaters
mliks <- sapply(model.sim, function(X){ X$mlik })
mliks <- mliks + sapply(b.sim, prior.beta)

probs <- exp(mliks - min(mliks))
probs <- probs / sum(probs)


save(file = "INLA-MHbiv.RData", list = ls())

#BMA models
library(INLABMA)

#bma.model <- INLABMA(model.sim[-(1:n.burnin)], 1:(n.sim-n.burnin))

models <- lapply(model.sim, function(X){X$model})
ws <- rep(1/length(models), length(models))

listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models, ws, X)
    })


#Bandwidth for densities 
bandw <- 0.08

pdf(file = "INLA-MHbiv2.pdf", width = 7.5, height = 7.5)
par(plt = c(.15, .95, .15, .9))
par(mfcol = c(2,3))
#Contour plots: INLA+MCMC
#plot(b.sim[idx,], xlab = "", ylab = "", main = "Markov Chain")
contour(z.inla, lty = 1, xlab =  expression(beta[1]), ylab =  expression(beta[2]))
points(c(2, -2), col = "black", pch = 4, cex = 2, lwd = 3)

legend("topright", lty = 1:2, legend = c("INLA+MCMC", "MCMC"), 
  bty = "n", cex = .8)
legend("right", pch = 4, legend = "Value", bty = "n", cex = 0.8)

#Contour plot: MCMC
contour(z.mcmc, lty = 2, xlab =  expression(beta[1]), ylab =  expression(beta[2]))

points(c(2, -2), col = "black", pch = 4, cex = 2, lwd = 3)

legend("topright", lty = 1:2, legend = c("INLA+MCMC", "MCMC"), 
  bty = "n", cex = .8)
legend("right", pch = 4, legend = "Value", bty = "n", cex = 0.8)

#Density beta1
plot(density(b.sim[, 1], bw = bandw), main = expression(beta[1]),
  xlab = "")
lines(m1$marginals.fixed[[2]], lty = 2)
#JAGS
lines(density(jm2.samp$beta1[1,,], bw = bandw), lty = 3)

legend("topright", lty = 1:3, legend = c("INLA+MCMC", "INLA", "MCMC"),
  bty = "n", cex = .8)

#Density beta2
plot(density(b.sim[, 2], bw = bandw), main = expression(beta[2]),
  xlab = "")
lines(m1$marginals.fixed[[3]], lty = 2)
#JAGS
lines(density(jm2.samp$beta2[1, , ], bw = bandw), lty = 3)

legend("topright", lty = 1:3, legend = c("INLA+MCMC", "INLA", "MCMC"),
  bty = "n", cex = .8)

#Fixed effects
plot(margeff[[1]][[1]], type = "l", xlab ="", ylab = "", main = expression(alpha))
lines(m1$marginals.fixed[[1]], lty = 2)
#JAGS
lines(density(jm2.samp$alpha[1,,], bw = 0.06), lty = 3)

legend("topright", lty = 1:3, legend = c("INLA+MCMC", "INLA", "MCMC"), bty = "n", cex = .8)

#Hyperparameters
plot(margeff[[2]][[1]], type = "l", xlab ="", ylab = "", main = expression(tau))
lines(m1$marginals.hyper[[1]], lty = 2)
#JAGS
lines(density(jm2.samp$prec[1,,], bw = 0.02), lty = 3)

legend("topright", lty = 1:3, legend = c("INLA+MCMC", "INLA", "MCMC"),
  bty = "n", cex = .8)


dev.off()
