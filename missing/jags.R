#Run mdoels with JAGS
library(rjags)

#Load data
#y <- 3 + 2*x+err
#err <- rnorm(n)
load("data_miss.RData")

d.mis <- as.list(d.mis)

d.mis$N <- length(d.mis$y)

d.mis$idx.mis <- which(is.na(d.mis$x))
d.mis$n.mis <- length(d.mis$idx.mis)

#Mean and precision for prior on the missing values
d.mis$mean.mis <- mean(d.mis$x, na.rm = TRUE)
d.mis$prec.mis <- 1/(4*var(d.mis$x, na.rm = TRUE))#2 times the s.d.

#MISSING VALUES MODEL
jm1 <- jags.model('modelmis.bug',
                   data = d.mis,
                   n.chains = 1,
                   n.adapt = 100)

update(jm1, 500)

jm1.samp <- jags.samples(jm1,
  c('alpha', 'beta', 'prec', 'mu', 'x'),
  n.iter = 100000, thin = 10)


print(jm1.samp)

save(file = "jags.RData", list = c("jm1.samp"))

