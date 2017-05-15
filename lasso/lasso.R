#
#Example on how to combine MCMC and R-INLA for Bayesian lasso
#

#Examples taken from http://www-bcf.usc.edu/~gareth/ISL/index.html, Ch. 6

#Load libraries
library(ISLR)
library(glmnet)

#
#Load data
#

data(Hitters)
summary(Hitters)

#Check NA's and fix
sum(is.na(Hitters$Salary))
Hitters <- na.omit(Hitters)

#
# The Lasso
#

#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)


#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid)
plot(lasso.mod)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1)
plot(cv.out)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])
mean((lasso.pred - y[test])^2)

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)

#Check estimated coefficients
lasso.coef
lasso.coef[lasso.coef != 0]

#Fitted values
lasso.fitted <- predict(out, s = bestlam, newx = x)



#Fit model with JAGS
library(rjags)

d.jags <- list(y = as.vector(y), x = x, p = 5, n = nrow(y),
  lambda = 1/0.073)
init.jags <- list(tau = 1, b = rep(0, 5))

lasso.m <- jags.model("lasso.bug", d.jags)
#Burn-in
jags.samples(lasso.m, c("b", "tau"), 500)
#Samples
smp.jags <- jags.samples(lasso.m, c("b", "tau"), 
  n.iter = 100000, thin = 10)
smp.jags

apply(smp.jags$b[,,1], 1, quantile)

lasso.coef


save(file = "lasso.RData", list = ls())


