##############################################
# Clean up and load relevant packages
rm(list=ls())
library(R2OpenBUGS)
library(coda)
library(mcmcplots)
library(dplyr)
library(ggplot2)
library(bayestestR)
library(arm)
library(tidyverse)
library(LearnBayes)
library(epiR)
library(data.table)
##############################################
#####################################################
# Dose Response Model
#####################################################
##############################################
# DATA PREPARATION
# Dose level
d <- c(0,62.5,125,250,500)
# Number of fetus
Nfetus <- c(282,225,290,261,141)
# Number of malformations
y <- c(67,34,193,250,141)
# Number of records
n <- 5
dyme1 <- data.frame(cbind(d,Nfetus,y))
##############################################
##############################################
# Run a logistic regression
dyme2 <- data.frame(list(malform=c(rep(1, 67),rep(0, 282-67),
                                   rep(1, 34),rep(0, 225-34),
                                   rep(1, 193),rep(0, 290-193),
                                   rep(1, 250),rep(0, 261-250),
                                   rep(1,141),rep(0, 141-141)),
                         dose=c(rep(0, 282),rep(62.5, 225), rep(125, 290), rep(250, 261), rep(500, 141))))
(glmfit <- glm(malform ~ dose, data = dyme2, family = binomial))
#alpha = -1.78190
#beta = 0.01823
##############################################
##############################################
#sampled priors
par(mfrow=c(2,2))
n_sample<-100000
prior_norm<-rnorm(n_sample, 0, 1000)
prior_norm_1<-rnorm(n_sample, 0, 1)
prior_t<-rt(n_sample, 4)
prior_logis<-rlogis(n_sample, 0, 1)
#plot pi density
plot(density(invlogit(prior_norm_1)), main = 'Normal distribution with variance 1')
plot(density(invlogit(prior_norm)), main = 'Normal distribution with variance 1000')
plot(density(invlogit(prior_t)), main = 't-distribution with 4 df')
plot(density(invlogit(prior_logis)), main = 'logistic 1 with scale 1')
par(mfrow=c(1,1))
##############################################
##############################################
# MCMC settings
ni <- 10000 #iterations
nb <- 5000 #burn-in
nc <- 2 #chains
##############################################
##############################################
# QUESTION NUMBER 1
# 1.1 - Write a bayesian model (likelihood and prior) using vague prior,
# 2 chains with 2 different sets of initials
model.data <- list("n"=n,'Nfetus' = Nfetus, 'y' = y, 'd'= d)
# Model : initial values
# since our model has alpha and beta, we specify sets of
# initials for these parameters, namely 2 sets of initials
# (any value).
model.inits2 <- list(alpha=0,beta=0) #equal to the mean of prior taken
model.inits1 <- list(alpha=glmfit$coef[1],beta=glmfit$coef[2])
model.inits = list(model.inits1, model.inits2)
# Specifying parameter names and OpenBUGS directory
parameters = c("alpha", "beta","ypred","bmd")
#parameters = c("alpha", "beta")
OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
# Model 1: Specification (Normal)
malformation <- function()
{
  for (i in 1:n)
  {
    y[i] ~ dbin(p[i],Nfetus[i])
    logit(p[i])<-alpha+beta*d[i]
  }
  # specifying priors - Remember that dnorm is specified as N(u,tau),
  # where tau is the precision of the normal
  # distribution = 1/sigma square. knowing this
  # we can specify a vague prior as a massive variance
  # which would be specifying an extremely small tau.
  alpha ~ dnorm(0, 0.0001)
  beta ~ dnorm(0, 0.0001)
  #Predict malformations at dose 100 for 240 exposed
  logit(pi2) <- alpha + beta * 100
  ypred ~ dbin(pi2,240)
  #Posterior distribution of the BMD
  P0 <- exp(alpha)/(1+exp(alpha))
  q.star <- (0.01*(1-P0)) + P0
  bmd <- (logit(q.star)-alpha)/beta
}
write.model(malformation, "model1.txt")
Sys.time()
# Model : running model
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits,
                  model.file = "model1.txt",
                  parameters=parameters,
                  n.chains = nc, n.iter = ni, n.burnin = nb,
                  codaPkg=T,
                  OpenBUGS.pgm=OpenBUGS.pgm)
Sys.time()
out1 <- read.bugs(model.out)
# summary statistics for posterior samples of the regression parameters
summary(out1)
HPDinterval(as.mcmc(as.matrix(out1)))
# History plot & posterior distributions
plot(out1, trace=TRUE, density = FALSE)
plot(out1, trace=FALSE, density = TRUE)
#Autocorrelation and running mean plots
acfplot(out1[,1:2], lag.max=50)
autocorr.plot(out1)
rmeanplot(out1,plot.title = "Running Mean Plots")
# Gelman-Rubin diagnostic
gelman.diag(out1)
gelman.plot(out1)
# Model 2: Specification (Logistic)
malformation <- function()
{
  for (i in 1:n)
  {
    y[i] ~ dbin(p[i],Nfetus[i])
    logit(p[i])<-alpha+beta*d[i]
  }
  # specifying priors
  alpha ~ dlogis(0, 0.1)
  beta ~ dlogis(0, 0.1)
  #Predict malformations at dose 100 for 240 exposed
  logit(pi2) <- alpha + beta * 100
  ypred ~ dbin(pi2,240)
  #Posterior distribution of the BMD
  P0 <- exp(alpha)/(1+exp(alpha))
  q.star <- (0.01*(1-P0)) + P0
  bmd <- (logit(q.star)-alpha)/beta
}
write.model(malformation, "model2.txt")
Sys.time()
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits,
                  model.file = "model2.txt",
                  parameters=parameters,
                  n.chains = nc, n.iter = ni, n.burnin = nb,
                  codaPkg=T,
                  OpenBUGS.pgm=OpenBUGS.pgm)
Sys.time()
out2 <- read.bugs(model.out)
# summary statistics for posterior samples of the regression parameters
summary(out2)
HPDinterval(as.mcmc(as.matrix(out2)))
# History plot & posterior distributions
plot(out2, trace=TRUE, density = FALSE)
plot(out2, trace=FALSE, density = TRUE)
#Autocorrelation and running mean plots
acfplot(out2[,1:2], lag.max=50)
autocorr.plot(out2)
rmeanplot(out2,plot.title = "Running Mean Plots")
# Gelman-Rubin diagnostic
gelman.diag(out2)
gelman.plot(out2)
##############################################
##############################################
# QUESTION NUMBER 2
# Model 3: Specification (t-dist)
malformation <- function()
{
  for (i in 1:n)
  {
    y[i] ~ dbin(p[i],Nfetus[i])
    logit(p[i])<-alpha+beta*d[i]
  }
  # specifying priors
  alpha ~ dt(0, 0.0001, 4)
  beta ~ dt(0, 0.0001, 4)
  #Predict malformations at dose 100 for 240 exposed
  logit(pi2) <- alpha + beta * 100
  ypred ~ dbin(pi2,240)
  #Posterior distribution of the BMD
  P0 <- exp(alpha)/(1+exp(alpha))
  q.star <- (0.01*(1-P0)) + P0
  bmd <- (logit(q.star)-alpha)/beta
}
write.model(malformation, "model3.txt")
Sys.time()
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits,
                  model.file = "model3.txt",
                  parameters=parameters,
                  n.chains = nc, n.iter = ni, n.burnin = nb,
                  codaPkg=T,
                  OpenBUGS.pgm=OpenBUGS.pgm)
Sys.time()
out3 <- read.bugs(model.out)
# summary statistics for posterior samples of the regression parameters
summary(out3)
HPDinterval(as.mcmc(as.matrix(out3)))
# History plot & posterior distributions
plot(out3, trace=TRUE, density = FALSE)
plot(out3, trace=FALSE, density = TRUE)
#Autocorrelation and running mean plots
acfplot(out3[,1:2], lag.max=50)
autocorr.plot(out3)
rmeanplot(out3,plot.title = "Running Mean Plots")
# Gelman-Rubin diagnostic
gelman.diag(out3)
gelman.plot(out3)
##############################################
##############################################
# QUESTION NUMBER 3
# See results from above (plots and summary statistics)
#Normal prior
plot(out1, trace=FALSE, density = TRUE)
# To assess the accuracy, we use the rule of thumb that our MC Error
#(Time-Series: Since it accounts for autocorrelation) is less than 5% of the standard deviation.
p1 <- summary(out1)[1]
p1 <- reduce(p1, bind_rows)
p1 <- as.data.frame(p1)
p1$`TSSE/SD` <-(p1$`Time-series SE`)/p1$`SD`
round(p1,3)
round(HPDinterval(as.mcmc(as.matrix(out1))),3)
#Logistic prior
plot(out2, trace=FALSE, density = TRUE)
# To assess the accuracy, we use the rule of thumb that our MC Error
#(Time-Series: Since it accounts for autocorrelation) is less than 5% of the standard deviation.
p2 <- summary(out2)[1]
p2 <- reduce(p2, bind_rows)
p2 <- as.data.frame(p2)
p2$`TSSE/SD` <-(p2$`Time-series SE`)/p2$`SD`
round(p2,3)
round(HPDinterval(as.mcmc(as.matrix(out2))),3)
#t-distribution prior
plot(out3, trace=FALSE, density = TRUE)
# To assess the accuracy, we use the rule of thumb that our MC Error
#(Time-Series: Since it accounts for autocorrelation) is less than 5% of the standard deviation.
p3 <- summary(out3)[1]
p3 <- reduce(p3, bind_rows)
p3 <- as.data.frame(p3)
p3$`TSSE/SD` <-(p3$`Time-series SE`)/p3$`SD`
round(p3,3)
round(HPDinterval(as.mcmc(as.matrix(out3))),3)
##############################################
##############################################
# QUESTION NUMBER 4
# normal priors
summary(out1)
# Store the chains in a data frame
weight_chains1 <- data.frame(out1[[1]])
weight_chains2 <- data.frame(out1[[2]])
weight_chains <- rbind(weight_chains1,weight_chains2)
(alpha1 <- mean(weight_chains$alpha))
(beta1 <- mean(weight_chains$beta))
(ypred1 <- mean(weight_chains$ypred))
(bmd1 <- mean(weight_chains$bmd))
#logistic priors
summary(out2)
# Store the chains in a data frame
weight_chains1 <- data.frame(out2[[1]])
weight_chains2 <- data.frame(out2[[2]])
weight_chains <- rbind(weight_chains1,weight_chains2)
(alpha2 <- mean(weight_chains$alpha))
(beta2 <- mean(weight_chains$beta))
(ypred2 <- mean(weight_chains$ypred))
(bmd2 <- mean(weight_chains$bmd))
#t-distribution priors
summary(out3)
# Store the chains in a data frame
weight_chains1 <- data.frame(out3[[1]])
weight_chains2 <- data.frame(out3[[2]])
weight_chains <- rbind(weight_chains1,weight_chains2)
(alpha3 <- mean(weight_chains$alpha))
(beta3 <- mean(weight_chains$beta))
(ypred3 <- mean(weight_chains$ypred))
(bmd3 <- mean(weight_chains$bmd))
# compute probabilities
dyme1 <- dyme1 %>%
  mutate(prob1 = exp(alpha1 + beta1*d)/(1+exp(alpha1 + beta1*d)),
         prob2 = exp(alpha2 + beta2*d)/(1+exp(alpha2 + beta2*d)),
         prob3 = exp(alpha3 + beta3*d)/(1+exp(alpha3 + beta3*d)),
         prob = y/Nfetus
  )
#With logistic
colors <- c("Observed" = "red", "Posterior (Normal prior)" = "blue",
            "Posterior (Logistic prior)" = "green", "Posterior (t-distribution prior)" = "orange")
ggplot(dyme1, aes(x=d)) +
  geom_line(aes(y = prob, color = "Observed"), size=1.2) +
  geom_line(aes(y = prob1, color="Posterior (Normal prior)"), size=1.1, linetype="twodash") +
  geom_line(aes(y = prob2, color="Posterior (Logistic prior)"), size=1.1, linetype="twodash") +
  geom_line(aes(y = prob3, color="Posterior (t-distribution prior)"), size=1.1, linetype="twodash") +
  labs(x="Dose", y="Probability", color = "Legend") +
  scale_color_manual(values = colors) +
  ggtitle("Observed probabilities vs Posterior dose-response relationship")
colors <- c("Observed" = "red", "Posterior (Normal prior)" = "blue",
            "Posterior (t-distribution prior)" = "orange")
#Without logistic
ggplot(dyme1, aes(x=d)) +
  geom_line(aes(y = prob, color = "Observed"), size=1.2) +
  geom_line(aes(y = prob1, color="Posterior (Normal prior)"), size=1.1, linetype="dotted") +
  geom_line(aes(y = prob3, color="Posterior (t-distribution prior)"), size=1.1, linetype="dashed") +
  labs(x="Dose", y="Probability", color = "Legend") +
  scale_color_manual(values = colors) +
  ggtitle("Observed probabilities vs Posterior dose-response relationship")
##############################################
##############################################
# QUESTION NUMBER 5
# Using normal priors
round(bmd1,3)
round(HPDinterval(as.mcmc(as.matrix(out1))),3)
#Using logistic priors
round(bmd2,3)
round(HPDinterval(as.mcmc(as.matrix(out2))),3)
#Using t-distribution priors
round(bmd3,3)
round(HPDinterval(as.mcmc(as.matrix(out3))),3)
#Safe level of exposure: around 3.74 but to be safe, will use lower bound which is 3.250
##############################################
##############################################
# QUESTION NUMBER 6
# Using normal priors
round(ypred1,3)
round(HPDinterval(as.mcmc(as.matrix(out1))),3)
#Using logistic priors
round(ypred2,3)
round(HPDinterval(as.mcmc(as.matrix(out2))),3)
#Using t-distribution priors
round(ypred3,3)
round(HPDinterval(as.mcmc(as.matrix(out3))),3)
#122 malformations for 240 exposed fetus at dose 100
#Corresponds to 0.508 which looking at the plot from question 4 is more or less the same.
##############################################
#####################################################
# Animal Prevalence
#####################################################
##############################################
# Data Preparation
dat1<-data.frame(N=c(272,87,322,176,94,387,279,194,65,110,266,397,152,231))
dat1$Z<-c(17,15,71,17,9,23,78,59,47,34,43,57,29,17)
dat1$type<-c(1,1,1,1,1,1,0,0,0,0,0,0,0,0)
##############################################
##############################################
# QUESTION NUMBER 1
# 2.1 - Write a bayesian model (likelihood and prior) using vague prior,
# 2 runned chains with 2 different sets of initials
model.data <- list("n"=14,'N' = dat1$N, 'Z' = dat1$Z, 'type'=dat1$type)
# Model : initial values
# since our model has beta0 and beta1, we specify sets of
# initials for these parameters, namely 2 sets of initials
# (any value).
model.inits <- list(list(beta0=0,beta1=0), list(beta0=1,beta1=1))
# function(){list(beta0=c(0,0)),beta1=c(0,1)}
# Model : Specification
#
animal <- function(){
  for (i in 1:n) {
    logit(p[i])<-beta0+beta1*type[i]
    Z[i] ~ dbin(p[i],N[i])
  }
  # specifying priors - Remember that dnorm is specified as N(u,tau),
  # where tau is the precision of the normal
  # distribution = 1/sigma square. knowing this
  # we can specify a vague prior as a massive variance
  # which would be specifying an extremely small tau.
  beta0 ~ dnorm(0, 0.00001) #
  beta1 ~ dnorm(0, 0.00001)
}
write.model(animal, "animal.txt")
# file.show("animal.txt")
# Specifying parameter names and OpenBUGS directory
parameters = c("p","beta0", "beta1")
OpenBUGS.pgm="C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
# Model : running model
# specify model, data, number of parallel chains
model.out <- bugs(model.data, model.inits,
                  model.file = "animal.txt",
                  parameters=parameters,
                  n.chains = 2, n.iter = 10000, n.burnin = 5000,
                  codaPkg=T,
                  OpenBUGS.pgm=OpenBUGS.pgm)
##############################################
##############################################
# QUESTION NUMBER 2
# 2.2 - Run 5000 iterations then look at history plots and autocorrelation,
# do a Gelman and Rubin convergence diagnostic and check convergence
# Model : running model with 5000 iterations
#
model.out2 <- bugs(model.data, model.inits,
                   model.file = "animal.txt",
                   parameters=parameters,
                   n.chains = 2, n.iter = 5000, n.burnin = 2500,
                   codaPkg=T,
                   OpenBUGS.pgm=OpenBUGS.pgm)
# doing history plots and showing posterior distribution
out <- read.bugs(model.out2)
par(mfrow=c(1,2))
traceplot(out[,1:2])
# doing autocorrelation plots
par(mfrow=c(1,1))
acfplot(out[,1:2], lag.max=50)
# Gelman and rubin convergence diagnostic
combinedchains = mcmc.list(out[[1]], out[[2]])
par(mfrow=c(1,2))
gelman.plot(combinedchains[,1:2])
# simulations seem they have converged because have a long tail close to 1
##############################################
##############################################
# QUESTION NUMBER 3
# 2.3 - Generate summary statistics and kernel density plots for hte posterior
# samples of the regression parameters
# summary statistics for psoterior samples of the regression parameters
summary(out)
HPDinterval(as.mcmc(as.matrix(out)))
# plotting posterior distribution
par(mfrow=c(1,2))
densplot(out[,1:2])
# To assess the accuracy, we use the rule of thumb that our MC Error
#(Time-Series: Since it accounts for autocorrelation) is less than 5% of the standard deviation.
p <- summary(out)[1]
p <- reduce(p, bind_rows)
p <- as.data.frame(p)
p$`TSSE/SD` <-(p$`Time-series SE`)/p$`SD`
round(p, digits = 4)
# In our case the MCMC SE is about 1.6% of the posterior standard deviation
##############################################
##############################################
# QUESTION NUMBER 4
# 2.4 - Determine the apparent animal prevalence p based on this survey, and
# determine whether the 2 farms are different
# since this is a logit model, then we can determine difference
# with odd ratios:
# Dairy farm = type 0 -> exp(beta0+beta1*0)/1+exp(beta0+beta1*0)
b0 <- p$Mean[1]
b1 <- p$Mean[2]
pBeef <- p$Mean[4]
pDairy <- p$Mean[10]
(Dairy_farm<-exp(b0+(b1*0))/(1+exp(b0+b1*0)))
1/Dairy_farm # Dairy_farm is the same as pDairy
(1-pDairy)/pDairy # for every 3.654856 healthy cow there is 1 diseased cow
# Beef farm = type 1 -> exp(beta0+beta1*1)/1+exp(beta0+beta1*1)
(Beef_farm<-exp(b0+(b1*1))/(1+exp(b0+b1*1)))
1/Beef_farm # Beef_farm is the same as pBeef
(1-pBeef)/pBeef # for every 7.794579 healthy cow there is 1 diseased cow
(Aparent_prev <- dat1 %>% group_by(type) %>% summarise(p = Z/N) %>% summarise(mean(p)))
# Given that p explains the presence of certain disease, a higher value
# will mean greater presence of disease. Here we can see that the both
# farms are differente regarding the tested disease, for which the Dairy
# farm has worst number of cases (need to check)
#credibility intervals are not crossing, therefore we can conclude
#that they are different
(glmfit <- glm((dat1$Z/dat1$N) ~ type, data = dat1, family = binomial))
#
(Dairy_farm_glm<-exp(-0.9785+(-1.0343*0))/(1+exp(-0.9785-1.0343*0)))
1/Dairy_farm_glm
# Beef farm = type 1 -> exp(beta0+beta1*1)/1+exp(beta0+beta1*1)
(Beef_farm_glm<-exp(-0.9785+(-1.0343*1))/(1+exp(-0.9785-1.0343*1)))
1/Beef_farm_glm
##############################################
##############################################
# QUESTION NUMBER 5
# 2.5 - plot posterior distributions to compare prevalence in the two types of farms
chainz <- rbindlist(lapply(out[,c(4,10)], as.data.frame))
names(chainz) <- c("Beef", "Dairy")
chainz2 <- melt(chainz)
(p <- ggplot(chainz2, aes(x = value, color = variable)) +
    geom_density(size=1.05) + geom_abline(colour='white',size = 1.2) + theme_bw())
# We can see both density plots separately as
densplot(out[,c(4,10)])
##############################################
##############################################
# QUESTION NUMBER 6
# 2.6 - Implementing priors given most probable values of the diagnostic sensitivity
# and specificity: 0.85 (95IC: 0.82-0.90), 0.95 (95IC: 0.9-0.97).
#
# Computing beta priors (optimization for alpha and beta parameters)
Se.Low.int=list(p=.025, x=0.82) # 2.5% quantile must be 0.82
Se.Up.int=list(p=.975, x=0.90) # 97.5% quantile must be 0.90
(Se.shape<-beta.select(Se.Low.int, Se.Up.int))
(Se.alpha0<-Se.shape[1])
(Se.beta0<-Se.shape[2])
(Se.mode<-(Se.alpha0-1) / (Se.alpha0+Se.beta0-2))
Sp.Low.int=list(p=.025, x=0.90) # 2.5% quantile must be 0.90
Sp.Up.int=list(p=.975, x=0.97) # 97.5% quantile must be 0.97
(Sp.shape<-beta.select(Sp.Low.int, Sp.Up.int))
(Sp.alpha0<-Sp.shape[1])
(Sp.beta0<-Sp.shape[2])
(Sp.mode<-(Sp.alpha0-1) / (Sp.alpha0+Sp.beta0-2))
# here we see that our mode is 86 and not 85 as should.
# So we tried "epi" package to compute the best possible
# parameters alpha and beta for our sensitivity prior.
(Se.shapes<-epi.betabuster(0.85, 0.95, greaterthan=T, x=0.825, conf.level = 0.95,
                           max.shape1 = 1000,step = 0.001))
#so we replace our alpha and beta for our optimal values
(Se.alpha0<-as.numeric(Se.shapes[1]))
(Se.beta0<-as.numeric(Se.shapes[2]))
#seeing the Se mode
(mode.o<-(Se.alpha0-1) / (Se.alpha0+Se.beta0-2))
#seeing the Se confidence interval
qbeta(p=c(0.025, 0.975), shape1=Se.alpha0, shape2=Se.beta0)
#Now we see that our prior fits really well our needs and we
# compute the parameters for the Specificity prior (Sp)
(Sp.shapes<-epi.betabuster(0.95, 0.95, greaterthan=T, x=0.911, conf.level = 0.95,
                           max.shape1 = 1000,step = 0.001))
(Sp.alpha0<-as.numeric(Sp.shapes[1]))
(Sp.beta0<-as.numeric(Sp.shapes[2]))
#seeing the Sp mode
(mode.o<-(Sp.alpha0-1) / (Sp.alpha0+Sp.beta0-2))
#seeing the Sp confidence interval
qbeta(p=c(0.025, 0.975), shape1=Sp.alpha0, shape2=Sp.beta0)
##############################################
##############################################
# QUESTION NUMBER 7
# 2.7 - Implementing priors given most probable values of the diagnostic sensitivity
# and specificity: 0.85 (95IC: 0.82-0.90), 0.93 (95IC: 0.9-0.97).
(tp1<-(0.1137+0.95-1)/(0.85+0.95+1))
(tp2<-(0.2148+0.95-1)/(0.85+0.95+1))
##############################################