## Example: Classification with INLA within MCMC

library(MASS)#For 'geyser' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
INLA:::inla.dynload.workaround()

library(INLABMA)
library(parallel)
options(mc.cores = 2)


#Get data
yy <- faithful$eruptions

#Number of data
n <- length(yy)

#PLot density
pdf(file = "geyser.pdf", width = 10, height = 5)
  par(mfrow = c(1, 2))

  plot(faithful$eruptions, faithful$waiting,
    xlab = "Eruption time", ylab = "Waiting time", pch = 19)
  plot(density(yy), xlab = "Eruption time", main = "")
dev.off()

## Initial split

#Number of groups
n.grp <- 2

# Using 3 as a cut-off point
#grp <- 1 + (yy >=3)

# 50-50 random assignment
#grp <- sample(rep(1:n.grp, length.out = n))

#1/3, 2/3's at random
#grp <- sample(1:2, n, rep = TRUE, prob = c(1/3, 2/3))
#Grouping using values
#grp <- 1 + as.integer(yy > 3)

#1/3, 2/3's in increasing order
grp <- rep(2, n)
grp[order(yy)[1:floor(n/3)]] <- 1


#The sampled values takes two parameters:
#z: the sampled index
#m: The fitted model

fit.inla <- function(yy, grp) {
  return(grp$m)
}

fit.inla.internal <- function(yy, grp) {

  #Data in two-column format
  y <- matrix(NA, ncol = n.grp, nrow = n)
  for(i in 1:n.grp) {
    idx <- which(grp == i)
    y[idx, i] <- yy[idx]
  }

  #X stores the intercept in the model
  x <- y
  x[!is.na(x)] <- 1
  
  d <- list(y = y, x = x)

  #Initial group
  m1 <- inla(y ~ -1 + x, data = d,
    family = rep("gaussian", n.grp),
    control.fixed = list(mean = list(x1 = 2, x2 = 4.5), prec = 1)
  )


  res<- list(model = m1, mlik = m1$mlik[1, 1])
  return(res)
}


#Completely at random
#x --> y
dq.z <- function(x, y, log = TRUE) {
  res <- log(0.5) * length(x)
  if(log) {
    return(res)
  }
  else {
    return(exp(res))
  }
}

rq.z <- function(z) {
  return(sample(1:n.grp, length(z), rep = TRUE))
}

#Probabilities of belonging to each group
get.probs <- function(z) {
  probs <- rep(0, n.grp)
  tab <- table(z) #+ 1 #Dirichlet prior
  probs[as.integer(names(tab))] <- tab/sum(tab)
  return(probs)
}

#Using means of conditional marginals
#FIXME: We do not consider possble label switching here
dq.z <- function(x, y, log = TRUE) {
  m.aux <- x$m$model #fit.inla(yy, x)$model
  means <- m.aux$summary.fixed[, "mean"]
  precs <- m.aux$summary.hyperpar[, "mean"]

  ww <- get.probs(x$z)

  z.probs <- sapply(1:n, function (X) {
    aux <- ww * dnorm(yy[X], means, sqrt(1/precs))
    (aux/sum(aux))[y$z[X]]
  })

  if(log) {
    return(sum(log(z.probs)))
  } else {
    return(prod(z.probs))
  }
}

#FIXME: We do not consider possible label switching here
rq.z <- function(z) {
  m.aux <- z$m$model #fit.inla(yy, z)$model#model.cur #
  means <- m.aux$summary.fixed[, "mean"]
  precs <- m.aux$summary.hyperpar[, "mean"]

  probs <- get.probs(z$z)

  #Sample
  #ww <- NA
  #ww <- (n * probs + 1)/(n + n.grp)
  #ww <- as.vector(rmultinom(1, n, ww)/n)

  #Re-sample probs
  #probs <- rmultinom(1, length(z$z), probs)
  #probs <- probs/sum(probs)

  z.sim <- sapply(1:n, function (X) {
    aux <- probs * dnorm(yy[X], means, sqrt(1/precs))
    sample(1:n.grp, 1, prob = aux/sum(aux))
  })

  #Fit model
  z.model <- fit.inla.internal(yy, z.sim)

  #Ne value
  z.new <- list(z = z.sim, m = z.model)#, w = ww)

  #Remove later
  #print("----------------------")
  #print(means)
  #print(precs)
  #print(ww)
  #print(z$m$mlik)
  #print(z.model$mlik)
  #print(dq.z(z.new, z) - dq.z(z, z.new))

  return(z.new)
}


prior.z <- function(z, log = TRUE) {

  res <- log(0.5) * length(z$z)
  if(log) {
    return(res)
  }
  else {
    return(exp(res))
  }
}

grp.init <- list(z = grp, m = fit.inla.internal(yy, grp),
 w = as.vector((table(grp) + 1) /(n + n.grp))
)

grp.init <- rq.z(grp.init)

#Run simulations
inlamh.res <- INLAMH(faithful$eruptions, fit.inla, grp.init, rq.z, dq.z, 
  prior.z, n.sim = 10000, n.burnin = 500, n.thin = 10, verbose = TRUE)

zz <- do.call(rbind, lapply(inlamh.res$b.sim, function(X){X$z}))

zz.probs <- apply(zz, 2, get.probs)

#Display probabilities
plot(zz.probs[1,], yy)

save(file = "oldfaith.RData", list = ls())

