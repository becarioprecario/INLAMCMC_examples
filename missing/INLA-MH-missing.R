#Implementation of INLA within M-H for Bayesian Inference
#Sampling of missing covariates

#Load libraries
library(INLA)
INLA:::inla.dynload.workaround()

#INLAMH()
source("INLAMH_function.R")

#Simulate data
n <- 100
set.seed(123)
x <- runif(n)
err <- rnorm(n)
y <- 3 + 2*x+err

d <- data.frame(y = y, x = x)

#INLA models
m1 <- inla(y ~ 1+x, data = d)

#Create missing covariates
n.mis <- 9
set.seed(321)
idx.mis <- sample(1:100, n.mis)
d.mis <- d
d.mis$x[idx.mis] <- NA

save(file = "data_miss.RData", list = c("d", "d.mis"))



#Fit linear model with R-INLA with a fixed beta
fit.inla <- function(data, x.mis) {#, idx.mis) {

   data$x[idx.mis] <- x.mis

   res <- inla(y ~1+x, data = data)

   return(list(mlik = res$mlik[1,1], model = res))
}


#Proposal x -> y
#density
dq.beta <- function(x, y, sigma = sqrt(1/12), log =TRUE) {
	res <- dnorm(y, mean = x, sd = sigma, log = log)

	if(log) {
		return(sum(res))
	} else {
		return(prod(res))
	}


}
#random
rq.beta <- function(x, sigma = sqrt(1/12) ) {
	rnorm(length(x), mean = x, sd = sigma)
}

#Prior for beta
prior.beta <- function(x, mu = mean(d.mis$x, na.rm = TRUE), 
   sigma = 2*sd(d.mis$x, na.rm = TRUE), log = TRUE) {
   res <- dnorm(x, mean = mu, sd= sigma, log = log)

	if(log) {
		return(sum(res))
	} else {
		return(prod(res))
	}
}


idx.mis <- which(is.na(d.mis$x))
#Run simulations
inlamh.res <- INLAMH(d, fit.inla, rep(mean(d.mis$x, na.rm = TRUE), 9),
  rq.beta, dq.beta, prior.beta, n.sim = 10000, n.burnin = 500, n.thin = 10)

x.sim <- do.call(rbind, inlamh.res$b.sim)
model.sim <- inlamh.res$model.sim

save(file = "INLA-MH-missing.RData", list = ls())


#Load JAGS results
#
#RUn source("jags.R", echo = TRUE) to get jags.RData
#
library(rjags)
load("jags.RData")


pdf(file = "INLA-MH-missing.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

for(i in 1:n.mis) {
  plot(x.sim[,i], type = "l", main = idx.mis[i], ylim = c(-1, 2) )
  abline(h = d$x[idx.mis[i]], col = "red")
  #JAGS
  lines(jm1.samp$x[idx.mis[i],,], col ="blue")
}
dev.off()


#Density
pdf(file = "INLA-MH-missing-dens.pdf", width = 10, height = 10)
par(mfrow = c(3,3))

for(i in 1:n.mis) {
  plot(density(x.sim[,i], bw = 0.055), type = "l", 
    xlab = "",
    main = paste("Observation ",idx.mis[i]), ylim = c(-1, 2) )
  #JAGS
  lines(density(jm1.samp$x[idx.mis[i],,], bw = 0.055), lty = 2)
  #Value
  abline(v = d$x[idx.mis[i]], lty = 3)

  legend("topleft", lty = 1:3, legend = c("INLA+MCMC", "MCMC", "Value"), 
    cex = .75, bty = "n")
}

dev.off()


print(c("Acceptance rate:", mean(inlamh.res$acc.sim))) 


#Summary statistics
#print("Summary statistics of beta:")
#c(mean(b.sim[idx]), sd(b.sim[idx]), quantile(b.sim[idx], c(0.025, .5, .975)) )
#m1$summary.fixed[2,]

#We need to add the prior of the paramaters
mliks <- sapply(model.sim, function(X){ X$mlik })
mliks <- mliks + sapply(b.sim[idx], prior.beta)

probs <- exp(mliks-min(mliks))
probs <- probs/sum(probs)


save(file = "INLA-MH-missing.RData", list = ls())
