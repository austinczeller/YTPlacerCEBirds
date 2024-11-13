#simulations with lots of zeros
# https://aosmith.rbind.io/2019/03/06/lots-of-zeros/
  
  ##negative binomial, generalised poisson, zero inflated models
  
library(ggplot2) # v. 3.1.0
library(HMMpa) # a function that draws random samples from the poisson distribution
library(MASS) # glm.nb()

#negative binomial----

set.seed(16)
dat = data.frame(Y = rnbinom(200, mu = 10, size = .05) ) #mean (lambda) of 10, theta = .05. 
#variance of the neg. binom distribution = lambda + (lambda^2/theta). 
#small theta should have many zeros and a few large counts

ggplot(dat, aes(x = Y) ) +
  geom_histogram(binwidth = 5)  +
  theme_bw(base_size = 18) +
  labs(y = "Frequency",
       title = "Negative binomial",
       subtitle = "mean = 10, theta = 0.05" ) +
  annotate(geom = "text",
           label = paste("Proportion 0:", mean(dat$Y == 0), 
                         "\nMax Count:", max(dat$Y) ),
           x = 150, y = 100, size = 8)

#generalised poisson (not poisson!) ----
#seems to have more extreme maximum counts 
#when lambda2 = 0, same as poisson

set.seed(16)
dat = data.frame(Y = rgenpois(200, lambda1 = 0.5, lambda2 = 0.95) )

ggplot(dat, aes(x = Y) ) +
  geom_histogram(binwidth = 5)  +
  theme_bw(base_size = 18) +
  labs(y = "Frequency",
       title = "Generalized Poisson",
       subtitle = "lambda1 = 0.5, lambda2 = 0.95") +
  annotate(geom = "text",
           label = paste("Proportion 0:", mean(dat$Y == 0), 
                         "\nMax Count:", max(dat$Y) ),
           x = 600, y = 100, size = 8)

##lots of zeros or excess zeros? ----
#negative binomial & generalised poisson can have lots of zeros. 
##can check whether negative binomial is sufficient my checking the number of expected zeros
#for the model vs the number of actual zeros. 

set.seed(16)
x = runif(200, 5, 10) # simulate explanatory variable
b0 = 1 # set value of intercept
b1 = 0.25 # set value of slope
means = exp(b0 + b1*x) # calculate true means
theta = 0.25 # true theta
y = rnbinom(200, mu = means, size = theta) #simulated count

fit1 = glm.nb(y ~ x) #fit a model

#check for excess
sum(y == 0) #observed number of zeros  out of 200 counts (from data)

##check expected number of zeros for the distribution
preds = predict(fit1, type = "response") # estimated means
esttheta = summary(fit1)$theta # estimated theta

prop0 = dnbinom(x = 0, mu = preds, size = esttheta ) ## probability of zero for each pred?
round( sum(prop0) ) ##expect 72 zeros based on modeled distribution
##this model is a good fit (not surprising, because used nbinom to simulate data)

##an example with excess zeros ----

fit2 = glm(y ~ x, family = poisson)
sum(y == 0)

round( sum( dpois(x = 0,
                  lambda = predict(fit2, type = "response") ) ) )
##wouldn't expect any zeros if the distribution was poisson. 
##could also check model for over dispersion

AIC(fit1, fit2) ## fit 1 much better!

##Simulate overdispersed neg binom ----
set.seed(16)
x = runif(200, 5, 10) # simulate explanatory variable
b0 = 1 # set value of intercept
b1 = 0.25 # set value of slope
means = exp(b0 + b1*x) # calculate true means
theta = 0.25 # true theta
p0 <- 0.7 ##prob of zero, doesn't depend on covariates
y1 <- rnbinom(200, mu = means, size = theta)
y <- ifelse(rbinom(200, 1, p0) == 1, y1, 0) #simulated count

##if
fit1 = glm.nb(y ~ x) #fit a model

#check for excess
sum(y == 0) #observed number of zeros  out of 200 counts (from data)
preds = predict(fit1, type = "response") # estimated means
esttheta = summary(fit1)$theta # estimated theta

prop0 = dnbinom(x = 0, mu = preds, size = esttheta ) ## probability of zero for each pred?
round( sum(prop0) ) ##160 expected vs 161 actual... so fits??
fit1 ## but intercept estimate is lower than should be

hist(y)
hist(y1)

##where zero inflation depends on x----
set.seed(16)
x = runif(200, 5, 10) # simulate explanatory variable
b0 = 1 # set value of intercept
b1 = 0.25 # set value of slope
means = exp(b0 + b1*x) # calculate true means
theta = 0.25 # true theta
p0 = .3
p1 = .3
p0 <- plogis(p0 + p1*x) ##prob of zero

y1 <- rnbinom(200, mu = means, size = theta)
y <- ifelse(rbinom(200, 1, p0) == 1, y1, 0) #simulated count

fit1 = glm.nb(y ~ x) #fit a model

#check for excess
sum(y == 0) #observed number of zeros  out of 200 counts (from data)
preds = predict(fit1, type = "response") # estimated means
esttheta = summary(fit1)$theta # estimated theta

prop0 = dnbinom(x = 0, mu = preds, size = esttheta ) ## probability of zero for each pred?
round( sum(prop0) ) ##83 expected vs 79 actual... so fits??
fit1 ## 

hist(y)
hist(y1)

