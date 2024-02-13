# part 1 of code of ulam/rethinking tutorial (part 2)
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-2/
# this has the code part that sets up the models and fits them
# this code might take a good while to run

## ---- packages --------
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('rethinking') #for model fitting
library('fs') #for file path

## ---- data --------
# adjust as necessary
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
#using dataset 3 for fitting
#also removing anything in the dataframe that's not used for fitting
#makes the ulam/Stan code more robust
fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            dose_adj = simdat[[3]]$dose_adj,
            time = simdat[[3]]$time)
#pulling out number of observations
Nind = length(unique(simdat$m3$id))


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

## ---- model-2 --------
#narrow-prior, full-pooling model
#2x(N+2)+1 parameters
m2 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1*dose_adj,
  beta <-  b0[id] + b1*dose_adj,
  a0[id] ~ dnorm(mu_a,  0.0001),
  b0[id] ~ dnorm(mu_b, 0.0001),
  mu_a ~ dnorm(2, 1),
  mu_b ~ dnorm(0.5, 1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0,1)
)

## ---- model-3 --------
#regularizing prior, partial-pooling model
m3 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1*dose_adj,
  beta <-  b0[id] + b1*dose_adj,
  a0[id] ~ dnorm(2,  1),
  b0[id] ~ dnorm(0.5, 1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0,1)
)

## ---- model-4 --------
#adaptive priors, partial-pooling model
#2x(N+2)+1 parameters
m4 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1*dose_adj,
  beta <-  b0[id] + b1*dose_adj,
  a0[id] ~ dnorm(mu_a,  sigma_a),
  b0[id] ~ dnorm(mu_b, sigma_b),
  mu_a ~ dnorm(2, 1),
  mu_b ~ dnorm(0.5, 1),
  sigma_a ~ cauchy(0, 1),
  sigma_b ~ cauchy(0, 1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0, 1)
)

## ---- model-2a --------
#full-pooling model, population-level parameters only
#2+2+1 parameters
m2a <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0 + a1*dose_adj,
  beta <-  b0 + b1*dose_adj,
  a0 ~ dnorm(2,  0.1),
  b0 ~ dnorm(0.5, 0.1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0,1)
)

## ---- model-4a --------
#adaptive priors, partial-pooling model
#non-centered
m4a <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  #rewritten to non-centered
  alpha <-  mu_a + az[id]*sigma_a + a1*dose_adj,
  beta  <-  mu_b + bz[id]*sigma_b + b1*dose_adj,
  #rewritten to non-centered
  az[id] ~ dnorm(0, 1),
  bz[id] ~ dnorm(0, 1),
  mu_a ~ dnorm(2, 1),
  mu_b ~ dnorm(0.5, 1),
  sigma_a ~ cauchy(0, 1),
  sigma_b ~ cauchy(0, 1),
  a1 ~ dnorm(0.3, 1),
  b1 ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0, 1)

  )

## ---- model-5 --------
#no dose effect
#separate intercept for each individual/id
#2xN+1 parameters
m5 <- alist(
  # distribution of outcome
  outcome ~ dnorm(mu, sigma),

  # main equation for time-series trajectory
  mu <- exp(alpha)*log(time) - exp(beta)*time,

  #equations for alpha and beta
  alpha <-  a0[id],
  beta <-  b0[id],

  #priors
  a0[id] ~ dnorm(2,  10),
  b0[id] ~ dnorm(0.5, 10),

  sigma ~ cauchy(0,1)
)




## ---- startvalues --------
## Setting starting values
#starting values for model 1
startm1 = list(a0 = rep(2,Nind), b0 = rep(0.5,Nind), a1 = 0.3 , b1 = -0.3, sigma = 1)
#starting values for model 2
startm2 = list(a0 = rep(2,Nind), b0 = rep(0.5,Nind), mu_a = 2, mu_b = 1, a1 = 0.3 , b1 = -0.3, sigma = 1)
#starting values for model 3
startm3 = startm1
#starting values for models 4 and 4a
startm4 = list(mu_a = 2, sigma_a = 1, mu_b = 1, sigma_b = 1, a1 = 0.3 , b1 = -0.3, sigma = 1)
startm4a = startm4
#starting values for model 2a
startm2a = list(a0 = 2, b0 = 0.5, a1 = 0.3, b1 = -0.3, sigma = 1)
#starting values for model 5
startm5 = list(a0 = rep(2,Nind), b0 = rep(0.5,Nind), sigma = 1)

#put different starting values in list
#need to be in same order as models below
startlist = list(startm1,startm2,startm3,startm4,startm2a,startm4,startm5)



## ---- fittingsetup --------
#general settings for fitting
#you might want to adjust based on your computer
warmup = 3000 
iter = 2000
max_td = 15 #tree depth
adapt_delta = 0.999
chains = 5
cores  = chains
seed = 4321
# for quick testing, use the settings below
# results won't make much sense, but can make sure the code runs
#warmup = 600 #for testing
#iter = warmup + floor(warmup/2)
#max_td = 10 #tree depth
#adapt_delta = 0.99


#stick all models into a list
modellist = list(m1=m1,m2=m2,m3=m3,m4=m4,m2a=m2a,m4a=m4a,m5=m5)
# set up a list in which we'll store our results
fl = vector(mode = "list", length = length(modellist))


#setting for parameter constraints
constraints = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")


## ---- modelfitting --------
# fitting models
#loop over all models and fit them using ulam
for (n in 1:length(modellist))
{

  cat('************** \n')
  cat('starting model', names(modellist[n]), '\n')

  tstart=proc.time(); #capture current time

  #run model fit
  fit <- ulam(flist = modellist[[n]],
                          data = fitdat,
                          start=startlist[[n]],
                          constraints=constraints,
                          log_lik=TRUE, cmdstan=TRUE,
                          control=list(adapt_delta=adapt_delta,
                                       max_treedepth = max_td),
                          chains=chains, cores = cores,
                          warmup = warmup, iter = iter,
                          seed = seed
  )# end ulam

  # save fit object to list
  fl[[n]]$fit <- fit
  
  #capture time taken for fit
  tdiff=proc.time()-tstart;
  runtime_minutes=tdiff[[3]]/60;

  cat('model fit took this many minutes:', runtime_minutes, '\n')
  cat('************** \n')

  #add some more things to the fit object
  fl[[n]]$runtime = runtime_minutes
  fl[[n]]$model = names(modellist)[n]

} #end fitting of all models

# saving the list of results so we can use them later
# the file is too large for standard Git/GitHub
# Git Large File Storage should be able to handle it
# I'm using a simple hack so I don't have to set up Git LFS
# I am saving these large file to a folder that is synced with Dropbox
# adjust accordingly for your setup
#filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
saveRDS(fl,filepath)


