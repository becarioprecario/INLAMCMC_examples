#Run mdoels with JAGS
library(rjags)

#Load data
#y <- 3 + 2*x+err
#err <- rnorm(n)
load("data.RData")

d <- as.list(d)

d$N <- length(d$y)

#1 COVARIATE DATA
jm1 <- jags.model('model.bug',
                   data = d,
                   n.chains = 1,
                   n.adapt = 100)

update(jm1, 500)

jm1.samp <- jags.samples(jm1,
             c('alpha', 'beta', 'prec', 'mu'),
             n.iter = 100000, thin = 10)

print(jm1.samp)


#TWO COVARIATES DATA
#y <- 3 + 2*x1-2*x2+err
#err <- rnorm(n)
load("databiv.RData")

d <- as.list(d)
d$N <- length(d$y)

jm2 <- jags.model('modelbiv.bug',
                   data = d,
                   n.chains = 1,
                   n.adapt = 100)

update(jm2, 500)

jm2.samp <- jags.samples(jm2,
             c('alpha', 'beta1', 'beta2', 'prec', 'mu'),
             n.iter = 100000, thin = 10)

print(jm2.samp)


save(file = "jags.RData", list = c("jm1.samp", "jm2.samp"))


