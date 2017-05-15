## Example: Classification with INLA within MCMC

library(MASS)#For 'geyser' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
library(INLABMA)
library(parallel)
options(mc.cores = 2)
library(rjags)


#Get data
yy <- faithful$eruptions

#Number of data
n <- length(yy)

n.grp <- 2

#Create jags data
jags.data <- list(y = yy, N = length(yy),
  K = n.grp, alpha = rep(1, n.grp),
  mean.grp = c(2, 4.5), prec.grp = c(1, 1),
  gamma.a = 1, gamma.b = 5e-05)


#Initial grouping
grp <- rep(2, n)
grp[order(yy)[1:floor(n/3)]] <- 1

jags.inits <- list(Z = grp, #Z = sample(1:n.grp, jags.data$N, TRUE),
  mu.orig = seq(min(jags.data$y), max(jags.data$y), length.out = n.grp)
  #, sigma = rep(1, n.grp)
)

# This model uses Gamma priors on the precisions and
# the results match thoser obtained with INLA
m.jags <- jags.model("mixture.bug", jags.data, jags.inits)
update(m.jags, 200)
res.jags <- jags.samples(m.jags, c("mu", "prec", "Z", "w"),
  n.iter = 10000, thin = 10 )

res.jags

save(file = paste0("jags-gaussian-ngrp-", n.grp, ".RData"), list = ls())
