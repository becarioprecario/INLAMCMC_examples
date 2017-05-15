#Implementation of INLA within M-H for Bayesian Inference
#  SPATIAL ERROR MODEL (see INLABMA::sem.inla)

#Load libraries
library(INLA)

INLA:::inla.dynload.workaround()


library(spdep)
library(INLABMA)

load("INLA-MH-speco.RData")

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
sembma<-INLABMA(seminla, rrho, 0, usenormal = TRUE)


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

rm(inlamh.res)
rm(model.sim)
rm(b.sim)
rm(models)
save(file = "summary-speco.RData", list = ls())


