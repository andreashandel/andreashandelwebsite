## Setup info
# Windows 11, R 4.3.2, 
#rstan_2.32.5       StanHeaders_2.32.5 digest_0.6.34      
#fs_1.6.3           rethinking_2.40    posterior_1.5.0   
#cmdstanr_0.7.1     ggplot2_3.4.4      dplyr_1.1.4     
# cmdstan-2.34.0


## ---- packages --------
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('rethinking') #for model fitting


## ---- settings --------
## General settings
set.seed(123) #for reproducibility
# days at which we assume outcome is measured
timevec <- c(0.1,1,3,5,7,10,14,21,28,35,42)

#different number of individuals per dose to make it clearer which is which
#also, that's the structure of the data which motivated the tutorial
Nlow = 7; Nmed = 8; Nhigh = 9; filename = "simdat.Rds"
#if you want to explore how model fitting changes if you increase sample size
#turn on this line of code
#this is used in part 4 of the tutorial
#Nlow = 70; Nmed = 80; Nhigh = 90; filename = "simdat_big.Rds"

Ntot = Nlow + Nmed + Nhigh; #total number of individuals

# Set values for dose
# since we only consider dose on a log scale
# we'll log transform right here and then always use it in those log units
high_dose = log(1000)
med_dose = log(100)
low_dose = log(10)
dosevec = c(rep(low_dose,Nlow),rep(med_dose,Nmed),rep(high_dose,Nhigh))
# we are also creating a version of the dose variable
# that consists of ordered categories instead of numeric values
# we'll use that mostly for plotting
dosevec_cat = ordered(c(rep("low", Nlow),
                        rep("medium",Nmed),
                        rep("high",Nhigh)),
                      levels=c("low","medium","high"))


## Setting parameter values

## ---- commonpars --------
sigma = 1
a1 = 0.1
b1 = -0.1

## ---- m1pars --------
m1_mua = 3
m1_mub = 1
m1_sigmaa = 1
m1_sigmab = 1
m1_a0 = rnorm(n=Ntot, m1_mua, m1_sigmaa)
m1_b0 = rnorm(n=Ntot, m1_mub, m1_sigmaa)
# saving main parameters
m1pars = c(sigma = sigma, a1 = a1, b1 = b1,
           a0_mu = m1_mua, b0_mu = m1_mub)

## ---- m1sims --------
m1_alpha = m1_a0 + a1*(dosevec - med_dose)
m1_beta = m1_b0 + b1*(dosevec - med_dose)
#doing matrix multiplication to get time-series for each individual
#for that to work, the timevec vector needs to be transposed
m1_mu =  exp(m1_alpha) %*% t(log(timevec)) - exp(m1_beta) %*% t(timevec)
# apply variation following a normal distribution to each data point
m1_y = rnorm(length(m1_mu),m1_mu, sigma)
# in a final step, we reorganize the data into a long data frame with
# columns id, time, dose, model,
# the deterministic mean mu, and the normally distributed outcome.
# We store dose in 3 versions, the original (log transformed one),
# the one that has the middle value subtracted, and a categorical one.
# Note that trick using sort to get time in the right order.
# Not a robust way of doing things, but works here
m1_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)),
                     mu = as.vector(m1_mu),
                     outcome = as.vector(m1_y),
                     model = "m1")


## ---- makeplots --------
p1 <- ggplot(m1_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,200)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)") +
  theme_minimal()


#using dataset 3 for fitting
#also removing anything in the dataframe that's not used for fitting
#makes the ulam/Stan code more robust
fitdat=list(id=m1_dat$id,
            outcome = m1_dat$outcome,
            dose_adj = m1_dat$dose_adj,
            time = m1_dat$time)
#pulling out number of observations
Ntot = length(unique(m1_dat$id))


## ---- model-1 --------
#wide-prior, no-pooling model
#separate intercept for each individual/id
#2x(N+1)+1 parameters
m1 <- alist(
  # distribution of outcome
  outcome ~ dnorm(mu, sigma),
  
  # main equation for time-series trajectory
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  
  #equations for alpha and beta
  alpha <-  a0[id] + a1*dose_adj,
  beta <-  b0[id] + b1*dose_adj,
  
  #priors
  a0[id] ~ dnorm(2,  10),
  b0[id] ~ dnorm(0.5, 10),
  
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0,1)
)

## ---- startvalues --------
## Setting starting values
#starting values for model 1
startm1 = list(a0 = rep(2,Ntot), b0 = rep(0.5,Ntot), a1 = 0.3 , b1 = -0.3, sigma = 1)

## ---- fittingsetup --------
#general settings for fitting
#you might want to adjust based on your computer
warmup = 6000 
iter = warmup + floor(warmup/2)
max_td = 18 #tree depth
adapt_delta = 0.9999
chains = 5
cores  = chains
seed = 4321
# for quick testing, use the settings below
# results won't make much sense, but can make sure the code runs
warmup = 600 #for testing
iter = warmup + floor(warmup/2)
max_td = 10 #tree depth
adapt_delta = 0.99


#setting for parameter constraints
constraints = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")


#run model fit
fit <- ulam(flist = m1,
            data = fitdat,
            start=startm1, #if this is turned on, the returned fit object is wonky and the commands below fail
            constraints=constraints,
            log_lik=TRUE, cmdstan=TRUE,
            control=list(adapt_delta=adapt_delta,
                         max_treedepth = max_td),
            chains=chains, cores = cores,
            warmup = warmup, iter = iter,
            seed = seed
)# end ulam

print(fit)
trankplot(fit)

