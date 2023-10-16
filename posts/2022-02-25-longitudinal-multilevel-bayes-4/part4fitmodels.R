# code for fitting portion of tutorial part 4
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-4/

## ---- packages --------
library('fs') #for file path
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('rethinking') #for model fitting


## ---- model-2a --------
#full-pooling model, population-level parameters only
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

## ---- model-4 --------
#adaptive priors, partial-pooling model
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


## ---- model-4a --------
#adaptive priors, partial-pooling model
#no middle dose subtraction
m4a <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1 * dose,
  beta <-  b0[id] + b1 * dose,
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


## ---- model-5 --------
# different way of enforcing positive parameters
m5 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id]*(1 + (a2-1) * dose_adj2),
  beta <-  b0[id]*(1 + (b2-1) *dose_adj2),
  a0[id] ~ dlnorm(mu_a,  sigma_a),
  b0[id] ~ dlnorm(mu_b, sigma_b),
  mu_a ~ dnorm(1, 0.5),
  mu_b ~ dnorm(0, 0.5),
  sigma_a ~ cauchy(0, 1),
  sigma_b ~ cauchy(0, 1),
  a2 ~ dlnorm(0, 1),
  b2 ~ dlnorm(0, 1),
  sigma ~ cauchy(0, 1)
)



## ---- model-6 --------
#model with dose treated as categorical
m6 <- alist(
  outcome ~ dnorm(mu, sigma),
  mu <- exp(alpha)*log(time) - exp(beta)*time,
  alpha <-  a0[id] + a1[dose_cat2],
  beta <-  b0[id] + b1[dose_cat2],
  a0[id] ~ dnorm(mu_a,  sigma_a),
  b0[id] ~ dnorm(mu_b, sigma_b),
  mu_a ~ dnorm(2, 1),
  mu_b ~ dnorm(0.5, 1),
  sigma_a ~ cauchy(0, 1),
  sigma_b ~ cauchy(0, 1),
  a1[dose_cat] ~ dnorm(0.3, 1),
  b1[dose_cat] ~ dnorm(-0.3, 1),
  sigma ~ cauchy(0, 1)
)




## ---- fitsettings --------
######################################
# Define general fit settings
######################################
#general settings for fitting
#you might want to adjust based on your computer
warmup = 6000
iter = warmup + floor(warmup/2)
max_td = 18 #tree depth
adapt_delta = 0.9999
chains = 5
cores  = chains
seed = 123




## ---- modelfitting-function --------
######################################
# Function to fit each model
######################################
# function to run fit so I don't need to keep repeating
# some parameters are set globally.
# not very clean code but good enough for here :)
fitfunction <- function(model, data, start, constraints)
{
  tstart=proc.time(); #capture current time

  fl$fit <- ulam(flist = model,
                 data = data,
                 start=start,
                 constraints=constraints,
                 log_lik=TRUE, cmdstan=TRUE,
                 control=list(adapt_delta=adapt_delta, max_treedepth = max_td),
                 chains=chains, cores = cores,
                 warmup = warmup, iter = iter
  )# end ulam statement

  tend=proc.time(); #capture current time
  tdiff=tend-tstart;
  runtime_minutes=tdiff[[3]]/60;

  #add some more things to the fit object
  fl$runtime = runtime_minutes
  fl$model = names(model)

  return(fl)
}




## ---- fittingsetup-dat2 --------
######################################
# Model fit to alternative data set
######################################
#stick all models into a list
modellist = list(m2a=m2a, m4=m4)
# set up a list in which we'll store our results
fl = vector(mode = "list", length = length(modellist))

#now fitting dataset 2 we produced in the first post
#also removing anything in the dataframe that's not used for fitting
#makes the ulam/Stan code more robust
simdat <- readRDS("simdat.Rds")
fitdat = list(id=simdat[[2]]$id,
            outcome = simdat[[2]]$outcome,
            dose_adj = simdat[[2]]$dose_adj,
            dose_cat = simdat[[3]]$dose_cat,
            time = simdat[[2]]$time)
#pulling out number of observations
Ntot = length(unique(fitdat$id))

## Setting starting values
#starting values for model 2
startm2a = list(a0 = 2, b0 = 0.5, a1 = 0.5 , b1 = -0.5, sigma = 1)
#starting values for model 4
startm4 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.5 , b1 = -0.5, sigma = 1)
#put different starting values in list
#need to be in same order as models below
startlist = list(startm2a,startm4)

# defining constraints on parameters
constm2a = list(sigma="lower=0")
constm4 = list(sigma="lower=0")
constraintlist = list(constm2a,constm4)



## ---- modelfitting-dat2 --------
fl <- NULL
for (n in 1:length(modellist))
{
  cat('************** \n')
  cat('starting model', names(modellist[n]), '\n')
  fl[[n]] <- fitfunction(model = modellist[[n]],
                             data = fitdat,
                             start = startlist[[n]],
                             constraints = constraintlist[[n]])
  cat('model fit took this many minutes:', fl[[n]]$runtime, '\n')
  cat('************** \n')
}

# saving the results so we can use them later
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_dat2", ext="Rds")
saveRDS(fl,filepath)



## ---- fittingsetup-big --------
######################################
# Model fit to bigger data set
######################################
#stick all models into a list
modellist = list(m4=m4)
# set up a list in which we'll store our results
fits = vector(mode = "list", length = length(modellist))

#fitting dataset with larger sample size
simdat <- readRDS("simdat_big.Rds")
fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            dose_adj = simdat[[3]]$dose_adj,
            dose_cat = simdat[[3]]$dose_cat,
            time = simdat[[3]]$time
            )

#starting values for model 4
startm4 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.5 , b1 = -0.5, sigma = 1)
startlist = list(startm4)

# defining constraints on parameters
constm4 = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")
constraintlist = list(constm4)


## ---- modelfitting-big --------
# fitting model
fl <- NULL
cat('************** \n')
cat('starting model', names(modellist[1]), '\n')
fl[[1]] <- fitfunction(model = modellist[[1]],
                         data = fitdat,
                         start = startlist[[1]],
                         constraints = constraintlist[[1]])
cat('model fit took this many minutes:', fl[[1]]$runtime, '\n')
cat('************** \n')

# saving the results so we can use them later
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_big", ext="Rds")
saveRDS(fl,filepath)




## ---- fittingsetup-altpos --------
#stick all models into a list
modellist = list(m4a=m4a, m5=m5)
# set up a list in which we'll store our results
fits = vector(mode = "list", length = length(modellist))

#fitting dataset 3 we produced in the earlier post
#also removing anything in the dataframe that's not used for fitting
#makes the ulam/Stan code more robust
#note that we need both dose_adj and an alternative that divides by max dose
#for model 5
simdat <- readRDS("simdat.Rds")
fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            dose = simdat[[3]]$dose,
            dose_adj2 = simdat[[3]]$dose/max(simdat[[3]]$dose),
            dose_cat = simdat[[3]]$dose_cat,
            time = simdat[[3]]$time
            )
#pulling out number of observations
Ntot = length(unique(fitdat$id))


#starting values for model 4a
startm4a = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.5 , b1 = -0.5, sigma = 1)
#starting values for model 5
startm5 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.2 , b1 = -0.2, sigma = 1)
#put different starting values in list
#need to be in same order as models below
startlist = list(startm4a, startm5)

# defining constraints on parameters
constm4a = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")
constm5 = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")
constraintlist = list(constm4a, constm5)


## ---- modelfitting-altpos --------
# fitting models
fl <- NULL
for (n in 1:length(modellist))
{
  cat('************** \n')
  cat('starting model', names(modellist[n]), '\n')
  fl[[n]] <- fitfunction(model = modellist[[n]],
                         data = fitdat,
                         start = startlist[[n]],
                         constraints = constraintlist[[n]])
  cat('model fit took this many minutes:', fl[[n]]$runtime, '\n')
  cat('************** \n')
}
# saving the results so we can use them later
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_altpos", ext="Rds")
saveRDS(fl,filepath)



## ---- fittingsetup-cat --------
#stick all models into a list
modellist = list(m4=m4, m6=m6)
# set up a list in which we'll store our results

#note that we use dose_cat for plotting
#ulam wants categories coded as integers, so we make dose_cat2
simdat <- readRDS("simdat.Rds")
fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            dose = simdat[[3]]$dose,
            dose_adj = simdat[[3]]$dose_adj,
            dose_cat = simdat[[3]]$dose_cat,
            dose_cat2 = as.numeric(simdat[[3]]$dose_cat),
            time = simdat[[3]]$time)
#pulling out number of observations
Ntot = length(unique(fitdat$id))

#starting values
startm4 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.5 , b1 = -0.5, sigma = 1)
startm6 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = rep(0.3,3) , b1 = rep(-0.3,3), sigma = 1)
startlist = list(startm4, startm6)

# defining constraints on parameters
constm4 = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")
constm6 = list(sigma="lower=0",sigma_a="lower=0",sigma_b="lower=0")
constraintlist = list(constm4, constm6)


## ---- modelfitting-cat --------
# fitting model
fl <- NULL
for (n in 1:length(modellist))
{
  cat('************** \n')
  cat('starting model', names(modellist[n]), '\n')
  fl[[n]] <- fitfunction(model = modellist[[n]],
                         data = fitdat,
                         start = startlist[[n]],
                         constraints = constraintlist[[n]])
  cat('model fit took this many minutes:', fl[[n]]$runtime, '\n')
  cat('************** \n')
}

# saving the results so we can use them later
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_cat", ext="Rds")
saveRDS(fl,filepath)
cat('all done')

