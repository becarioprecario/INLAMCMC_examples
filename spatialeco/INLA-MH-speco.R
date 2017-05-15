#Implementation of INLA within M-H for Bayesian Inference
#  SPATIAL ERROR MODEL (see INLABMA::sem.inla)

#Load libraries
library(INLA)
library(spdep)
library(INLABMA)

INLA:::inla.dynload.workaround()
#INLAMH()
source("INLAMH_function.R")


#Load data and compute adjacency matrix
data(columbus)
lw <- nb2listw(col.gal.nb, style="W")

#FIT model using ML
colsemml <- errorsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
  quiet=FALSE)


#Adjacency matrix
W <- as(as_dgRMatrix_listw(nb2listw(col.gal.nb)), "CsparseMatrix")
#Index for spatial random effects
columbus$idx<-1:nrow(columbus)


     
#Formula
form<- CRIME ~ INC + HOVAL
     
zero.variance = list(prec=list(initial = 25, fixed=TRUE))


#Fit model for a given rho
fit.inla <- function(data, rho) {
#, form = form, W = W, 
#  zero.variance = zero.variance) {
  res <- sem.inla(form, d = data, W = W, rho = rho,
    family = "gaussian", impacts = FALSE,
    control.family = list(hyper = zero.variance),
    #control.predictor = list(compute = TRUE),
    #control.compute = list(dic = TRUE, cpo = TRUE),
    #control.inla = list(print.joint.hyper = TRUE), 
    ##tolerance=1e-20, h=1e-6),
    verbose = FALSE)

   return(list(mlik = res$mlik[1,1], model = res))
}


#Proposal x -> y
#density
dq.rho <- function(x, y, sigma = .15, log =TRUE) {
	dnorm(y, mean = x, sd = sigma, log = log)
}
#random
rq.rho <- function(x, sigma = .15) {
	rnorm(1, mean = x, sd = sigma)
}

#Prior for beta#Uniform [-1.5, 1] using eigenvalues
prior.rho <- function(x, log = TRUE) {
  dunif(x, -1.5, 1, log = log)
}

#Data as 'd'
d <- columbus

#Run simulations
inlamh.res <- INLAMH(d, fit.inla, 0, rq.rho, dq.rho, prior.rho,
  n.sim = 10000, n.burnin = 500, n.thin = 10)

#Show results
b.sim <- do.call(rbind, inlamh.res$b.sim)
model.sim <- inlamh.res$model.sim

save(file = "INLA-MH-speco.RData", list = ls())

#Compute model with JAGS
library(SEjags)

if(!file.exists("SEjags.RData")) {
  sem.mcmc <- SEjags(CRIME ~ INC + HOVAL, data = columbus, W = listw2mat(lw),
    model = "sem", n.burnin = 500, n.iter = 10000, n.thin = 10)

  save(file = "SEjags.RData", list = c("sem.mcmc"))
} else {
  load("SEjags.RData")
}

#Using Bivand et al. (2014)
library(parallel)
options(mc.cores = 4)

#Define grid on rho
rrho<-seq(-1, .95, length.out=40)

seminla<-mclapply(rrho, function(rho){
     
                     sem.inla(form, d=columbus, W=W, rho=rho,
                             family = "gaussian", impacts=FALSE,
                             control.family = list(hyper = zero.variance),
                             control.predictor=list(compute=TRUE),
                             control.compute=list(dic=TRUE, cpo=TRUE),
                             control.inla=list(print.joint.hyper=TRUE, 
                                     tolerance=1e-20, h=1e-6),
                             verbose=FALSE
                     )
     
})
sembma <- INLABMA(seminla, rrho, 0, usenormal = TRUE)


pdf(file = "INLA-MH.pdf", width = 10, height = 5)
par(mfrow = c(1,2))
plot(b.sim, type = "l")
abline(h = colsemml$lambda, col = "red")
abline(h = sembma$rho$mean, col = "red", lty = 2)
legend("topleft", legend = c("INLA+MCMC", "Bivand et al.", "Max. Lik."),
  lty = c(1, 2, 1), col = c("black", "red", "red"), bty = "n")

#Density
plot(density(b.sim))
lines(density(sem.mcmc$lambda), lty = 2)
lines(sembma$rho$marginal, lty = 3)
abline(v = colsemml$lambda, lty = 4)

legend("topleft", legend = c("INLA+MCMC", "MCMC", "Bivand et al.", "Max. Lik."),
  lty = 1:4, bty = "n")
#MCMC

dev.off()


print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 

#Summary statistics
print("Summary statistics of rho:")
c(mean(b.sim), sd(b.sim),
  quantile(b.sim, c(0.025, .5, .975)) )
inla.zmarginal(sembma$rho$marginal)

#We need to add the prior of the paramaters
mliks <- sapply(model.sim, function(X){ X$mlik })
mliks <- mliks + sapply(b.sim, prior.rho)

probs <- exp(mliks-min(mliks))
probs <- probs/sum(probs)


#BMA models
library(INLABMA)

#bma.model <- INLABMA(model.sim[-(1:n.burnin)], 1:(n.sim-n.burnin))

models <- lapply(model.sim, function(X){X$model})
n.sim <- length(models)
ws <- rep(1/n.sim, n.sim)

listmarg <- c("marginals.fixed", "marginals.hyperpar")
margeff <- mclapply(listmarg, function(X) {
        INLABMA:::fitmargBMA2(models, ws, X)
    })

#Summary estimates of fixed effects
do.call(rbind, lapply(margeff[[1]], inla.zmarginal))
sembma$summary.fixed

#Summary estimates
do.call(rbind, lapply(margeff[[2]], inla.zmarginal))
sembma$summary.hyperpar
sembma.hyper <- do.call(rbind, lapply(seminla, function(X){X$summary.hyperpar[,1:2]}))
mean(sembma.hyper[, 1]) #Mean
sqrt(mean(sembma.hyper[, 2]^2)) #Standard deviation


#MCMC summary
#-->rho
print(c("RHO MCMC (mean, sd)", mean(sem.mcmc$lambda[1,,1]), 
  sd(sem.mcmc$lambda[1,,1])))
#--->betas
apply(sem.mcmc$b[,1,,1], 1, mean)
#[1] 58.8026027 -0.9171711 -0.2978022
apply(sem.mcmc$b[,1,,1], 1, sd)
#[1] 6.4571607 0.4061325 0.1017372
#--->1/sigma^2_u
mean(sem.mcmc$tau[,,1])
sd(sem.mcmc$tau[,,1])






save(file = "INLA-MH-speco.RData", list = ls())



