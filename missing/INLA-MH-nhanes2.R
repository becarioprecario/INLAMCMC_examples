#Implementation of INLA within M-H for Bayesian Inference
#Sampling of missing covariates

library(mice)
data(nhanes2)


#Load libraries
library(INLA)
#  INLA:::inla.dynload.workaround()
#INLAMH()
library(INLABMA)

#INLA models
m1 <- inla(chl ~ 1 + bmi + age, data = nhanes2)

#Generic variables for model fitting
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)

save(file = "nhanes2.RData", list = c("nhanes2", "d.mis"))


#Fit linear model with R-INLA with a fixed beta
fit.inla <- function(data, x.mis) { #, idx.mis) {

   data$bmi[idx.mis] <- x.mis

   res <- inla(chl ~ 1 + bmi + age, data = data)

   return(list(mlik = res$mlik[1,1], model = res))
}


#Proposal x -> y
#density
dq.beta <- function(x, y, sigma = sqrt(10), log =TRUE) {
	res <- dnorm(y, mean = x, sd = sigma, log = log)

	if(log) {
		return(sum(res))
	} else {
		return(prod(res))
	}


}
#random
rq.beta <- function(x, sigma = sqrt(10) ) {
	rnorm(length(x), mean = x, sd = sigma)
}

#Prior for beta
prior.beta <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE), 
   sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
   res <- dnorm(x, mean = mu, sd= sigma, log = log)

	if(log) {
		return(sum(res))
	} else {
		return(prod(res))
	}
}



#Run simulations
inlamh.res <- INLAMH(d.mis, fit.inla, rep(mean(d.mis$bmi, na.rm = TRUE), n.mis),
  rq.beta, dq.beta, prior.beta, 
  n.sim = 10000, n.burnin = 500, n.thin = 10)

#Show results
x.sim <- do.call(rbind, inlamh.res$b.sim)
model.sim <- inlamh.res$model.sim

save(file = "INLA-MH-nhanes2.RData", list = ls())


#Load JAGS results
#
#Run source("jags.R", echo = TRUE) to get jags.RData
#
library(rjags)
load("jags.RData")

#Summary jags
mean(jm1.samp$alpha); sd(jm1.samp$alpha)
mean(jm1.samp$beta); sd(jm1.samp$beta)
mean(jm1.samp$b.age2); sd(jm1.samp$b.age2)
mean(jm1.samp$b.age3); sd(jm1.samp$b.age3)
mean(jm1.samp$prec); sd(jm1.samp$prec)

#This some of the observations
#	x.sim <- x.sim[1:8646, ]; model.sim <- model.sim[1:8646]

pdf(file = "INLA-MH-missing.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

for(i in 1:n.mis) {
  plot(x.sim[,i], type = "l", main = idx.mis[i])
  #JAGS
  lines(jm1.samp$bmi[idx.mis[i],,], col ="blue")
}
dev.off()


#Density
pdf(file = "INLA-MH-missing-dens.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

for(i in 1:n.mis) {
  plot(density(x.sim[,i], bw = 1.75), type = "l", 
    main = paste0("Observation ", idx.mis[i]), lty = 1,
  xlab = "")
  #JAGS
  lines(density(jm1.samp$bmi[idx.mis[i],,], bw = 1.75), col ="black", lty = 2)

  legend("topleft", lty = 1:2, legend = c("INLA within MCMC", "MCMC"),
    col = c("black", "black"), cex = .8, bty = "n")
}

dev.off()


print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 


#Summary statistics
#print("Summary statistics of beta:")
#c(mean(b.sim[idx]), sd(b.sim[idx]), quantile(b.sim[idx], c(0.025, .5, .975)) )
#m1$summary.fixed[2,]

#We need to add the prior of the paramaters
mliks <- sapply(model.sim, function(X){ X$mlik})
mliks <- mliks + apply(x.sim, 1, prior.beta)

probs <- exp(mliks - min(mliks))
probs <- probs/sum(probs)

#BMA models
library(INLABMA)
l.models <- lapply(model.sim, function(X){ X$model})
#m.final <- INLABMA(l.models, length(l.models), 0)

#Fixed effects
n.sim <- length(l.models)
ws <- rep(1/n.sim, n.sim)
INLABMA:::fitmatrixBMA(l.models, ws, "summary.fixed")
#> INLABMA:::fitmatrixBMA(l.models, ws, "summary.fixed")
#                 mean        sd  0.025quant  0.5quant 0.975quant      mode
#(Intercept) 43.471216 56.853046 -63.2769647 41.297175 162.666847 37.351921
#bmi          4.864194  2.007569   0.6950043  4.924653   8.686571  5.032377
#age40-59    29.500952 17.396997  -6.9732350 30.253022  61.739292 31.721176
#age60-99    49.448805 20.980671   5.1430652 50.527989  87.752216 52.721400
#                     kld
#(Intercept) 4.211303e-11
#bmi         2.744941e-11
#age40-59    5.957370e-11
#age60-99    3.403575e-11

#Hyperparameters
INLABMA:::fitmatrixBMA(l.models, ws, "summary.hyperpar")
#> INLABMA:::fitmatrixBMA(l.models, ws, "summary.hyperpar")
#                                              mean          sd   0.025quant
#Precision for the Gaussian observations 0.00108211 0.002121901 0.0004104344
#                                           0.5quant  0.975quant         mode
#Precision for the Gaussian observations 0.001005999 0.002113539 0.0008709942

#Marginals of fitted values
marg.fixed <- INLABMA:::fitmargBMA2(l.models, ws, "marginals.fixed")


warning("Here constaint to the interval 0.005 using debug().")
marg.hyperpar <- INLABMA:::fitmargBMA2(l.models, ws, "marginals.hyperpar")


#Summary statistics using inla.zmarginal
do.call(rbind, lapply(c(marg.fixed, marg.hyperpar), inla.zmarginal))

#Diplay fitted
pdf(file = "post_param.pdf", width = 7.5, height = 5)

par(mfrow = c(2, 3))
par(plt = c(.1, .95, .15, .85))


plot(marg.fixed[[1]], type = "l", main = expression(beta[0]), xlab = "")
lines(density(jm1.samp$alpha), lty = 2)
legend("topleft", legend = c("INLAMCMC", "MCMC"), lty = 1:2, bty = "n",
  cex = .7)

plot(marg.fixed[[2]], type = "l", main = expression(beta[1]), xlab = "")
lines(density(jm1.samp$beta), lty = 2)
legend("topleft", legend = c("INLAMCMC", "MCMC"), lty = 1:2, bty = "n",
  cex = .7)

plot(marg.fixed[[3]], type = "l", main = expression(beta[2]), xlab = "")
lines(density(jm1.samp$b.age2), lty = 2)
legend("topleft", legend = c("INLAMCMC", "MCMC"), lty = 1:2, bty = "n",
  cex = .7)

plot(marg.fixed[[4]], type = "l", main = expression(beta[3]), xlab = "")
lines(density(jm1.samp$b.age3), lty = 2)
legend("topleft", legend = c("INLAMCMC", "MCMC"), lty = 1:2, bty = "n",
  cex = .7)

plot(marg.hyperpar[[1]], type = "l", main = expression(tau), 
  xlim = c(0, .004), xlab = "" )
lines(density(jm1.samp$prec), lty = 2)
legend("topright", legend = c("INLAMCMC", "MCMC"), lty = 1:2, bty = "n",
  cex = .7)
dev.off()


save(file = "INLA-MH-nhanes2.RData", list = ls())
