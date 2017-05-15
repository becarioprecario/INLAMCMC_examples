#Run mdoels with JAGS
library(rjags)

#Load data
#y <- 3 + 2*x+err
#err <- rnorm(n)
load("data.RData")

d.mis <- as.list(d.mis)

#Age as categorical
d.mis$agecat <- model.matrix (~ -1 + age, data = data.frame(age = d.mis$age))

d.mis$N <- length(d.mis$chl)

d.mis$idx.mis <- which(is.na(d.mis$bmi))
d.mis$n.mis <- length(d.mis$idx.mis)

#Mean and precision for prior on the missing values
d.mis$mean.mis <- mean(d.mis$bmi, na.rm = TRUE)
d.mis$prec.mis <- 1/(4*var(d.mis$bmi, na.rm = TRUE))#2 times the s.d.

#MISSING VALUES MODEL
jm1 <- jags.model('modelmis.bug',
                   data = d.mis,
                   n.chains = 1,
                   n.adapt = 100)

update(jm1, 500)

jm1.samp <- jags.samples(jm1,
  c('alpha', 'beta', 'b.age2', 'b.age3', 'prec', 'mu', 'bmi', 'chl'),
  n.iter = 100000, thin = 10)

print(jm1.samp)

save(file = "jags.RData", list = c("jm1.samp"))

