#
#Implementation of Bayesian lasso with R-INLA and MCMC
#

#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()

#INLAMH()
library(INLABMA)

library(ISLR)
library(smoothmest)#Laplace distribution
library(mvtnorm)

#Fit linear model with R-INLA with a fixed beta
#data: list with 'y'y and 'x' (matrix of covariates, not including intercept)
fit.inla <- function(data, b) {

  data$oset <- data$x %*% matrix(b, ncol = 1)

  res <- inla(y ~ -1 + offset(oset), data = data)#,
#    control.inla = list(strategy = "laplace", int.strategy = "ccd",
#      dz = .1, restart = 10))

  res <- inla.rerun(res) #Double-check to get th right mlik
#  res <- inla.rerun(res) #Double-check to get th right mlik

  return(list(mlik = res$mlik[1,1], model = res))

  #Do NOT return the full model
#  return(list(mlik = res$mlik[1,1], model = NA))
}

#Prior for beta (each one is an independent Laplace distribution)
prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))

  if(!log) { res <- exp(res) }

  return(res)
}

#Conditional marg-likelihood (?)
cond.mlik <- function(b, data) {
   res <- - as.vector(fit.inla(data = data, b = b)$mlik + prior.beta (b, 0, 1/.08))

   print(c(res, b))

   return(res)
}

#Optim using a null starting point
#optim(rep(0, 5), cond.mlik, method = "BFGS", data = d)

#Optim using the solutio of the lasso
#optim(c(0, 1.818961, 0, 0, 4.052350), cond.mlik, method = "SANN", data = d)

#Load data
data(Hitters)
Hitters <- na.omit(Hitters)

#Data
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing: 1, 3, 4 should have a zero coef.
#x <- scale(x)#glmnet scales 'x' by default
y <- Hitters$Salary
#y <- y/sd(y) #glmnet sets 'y' yo have variance 1 (see manual page)

#Scale x and y
y <- scale(y)
x <- scale(x)
d <- list(y = y, x = x)

#Test code 
n.beta <- ncol(d$x)
summary( fit.inla(d, b = rep(0, n.beta))$model )

#LS estimate
x1 <-cbind(1,x)
ML.betas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
ML.betas
#            [,1]
#      164.097257
#AtBat  -1.802350
#Hits    6.308997
#HmRun  -3.431549
#Runs    2.038504
#RBI     6.745761
#b.sim[1, ] <- as.vector(ML.betas)[-1]

#lasso estimate
#b.sim[1, ] <- c(0, 1.782466e-01, 0, 0, 2.287608e-01)

#stdev.samp <- 1 * sqrt(diag(solve(t(x)%*%x))) * 
#  sd(y - x%*%solve(t(x)%*%x)%*%t(x)%*%y)
stdev.samp <- .25 * solve(t(x)%*%x)

#Proposal x -> y
#density
#dq.beta <- function(x, y, sigma = .01, log =TRUE) {
dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
        #sum(dnorm(y, mean = x, sd = sigma, log = log))
	dmvnorm(y, mean = x, sigma = sigma, log = log)
}
#random
#rq.beta <- function(x, sigma = .01) {
rq.beta <- function(x, sigma = stdev.samp) {
        #rnorm(length(x), mean = x, sd = sigma)
	as.vector(rmvnorm(1, mean = x, sigma = sigma))
}



#Run simulations
inlamh.res <- INLAMH(d, fit.inla, rep(0, n.beta), rq.beta, dq.beta, prior.beta,
  n.sim = 10000, n.burnin = 500, n.thin = 10, verbose = TRUE)

#Show results
b.sim <- do.call(rbind, inlamh.res$b.sim)
model.sim <- inlamh.res$model.sim

save(file = "INLA-lasso.RData", list = ls())
