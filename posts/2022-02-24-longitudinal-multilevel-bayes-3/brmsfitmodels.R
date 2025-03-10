# fitting portion of code of brms tutorial (part 3)
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-3/


## ---- packages --------
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('brms') # for model fitting
library('posterior') #for post-processing
library('fs') #for file path


## ---- data --------
simdat <- readRDS("simdat.Rds")
#pulling out number of observations
Ntot = length(unique(simdat$m3$id))

#fitting dataset 3
#we need to make sure the id is coded as a factor variable
#also removing anything in the dataframe that's not used for fitting
#makes the Stan code more robust
fitdat=list(id = as.factor(simdat[[3]]$id),
            outcome = simdat[[3]]$outcome,
            dose_adj = simdat[[3]]$dose_adj,
            time = simdat[[3]]$time)


## ---- model-1 --------
#no-pooling model
#separate intercept for each individual/id
#2x(N+1)+1 parameters
m1eqs <- bf(  #main equation for time-series trajectory
          outcome ~  exp(alpha)*log(time) - exp(beta)*time,
          #equations for alpha and beta
          alpha ~ 0 + id + dose_adj,
          beta  ~ 0 + id + dose_adj,
          nl = TRUE)

m1priors <- c(#assign priors to all coefficients related to both id and dose_adj for alpha and beta
              prior(normal(2, 10),  class = "b",  nlpar = "alpha"),
              prior(normal(0.5, 10),  class = "b",  nlpar = "beta"),
              #change the dose_adj priors to something different than the id priors
              prior(normal(0.3, 1),   class = "b",  nlpar = "alpha", coef = "dose_adj"),
              prior(normal(-0.3, 1),  class = "b",  nlpar = "beta", coef = "dose_adj"),
              prior(cauchy(0,1), class = "sigma") )



## ---- model-2a --------
#full-pooling model
#2+2+1 parameters
m2aeqs <- bf(  #main equation for time-series trajectory
  outcome ~ exp(alpha)*log(time) - exp(beta)*time,
  #equations for alpha and beta
  alpha ~ 1 + dose_adj,
  beta  ~  1 + dose_adj,
  nl = TRUE)

m2apriors <- c(prior(normal(2, 2),  class = "b",  nlpar = "alpha", coef = "Intercept"),
              prior(normal(0.5, 2),  class = "b",  nlpar = "beta", coef = "Intercept"),
              prior(normal(0.3, 1),   class = "b",  nlpar = "alpha", coef = "dose_adj"),
              prior(normal(-0.3, 1),  class = "b",  nlpar = "beta", coef = "dose_adj"),
              prior(cauchy(0,1), class = "sigma")  )


## ---- model-3 --------
#same as model 1 but regularizing priors
m3eqs <- m1eqs

m3priors <- c(#assign priors to all coefficients related to id and dose_adj for alpha and beta
  prior(normal(2, 1),  class = "b",  nlpar = "alpha"),
  prior(normal(0.5, 1),  class = "b",  nlpar = "beta"),
  #change the dose_adj priors to something different than the id priors
  prior(normal(0.3, 1),   class = "b",  nlpar = "alpha", coef = "dose_adj"),
  prior(normal(-0.3, 1),  class = "b",  nlpar = "beta", coef = "dose_adj"),
  prior(cauchy(0,1), class = "sigma") )


## ---- model-4 --------
#adaptive prior, partial-pooling model
m4eqs <- bf(  #main equation for time-series trajectory
  outcome ~ exp(alpha)*log(time) - exp(beta)*time,
  #equations for alpha and beta
  alpha ~  (1|id) + dose_adj,
  beta  ~  (1|id) + dose_adj,
  nl = TRUE)

m4priors <- c(prior(normal(2, 1),  class = "b",  nlpar = "alpha", coef = "Intercept"),
              prior(normal(0.5, 1),  class = "b",  nlpar = "beta", coef = "Intercept"),
              prior(normal(0.3, 1),   class = "b",  nlpar = "alpha", coef = "dose_adj"),
              prior(normal(-0.3, 1),  class = "b",  nlpar = "beta", coef = "dose_adj"),
              prior(cauchy(0,1), class = "sd", nlpar = "alpha"),
              prior(cauchy(0,1), class = "sd", nlpar = "beta"),
              prior(cauchy(0,1), class = "sigma")  )



## ---- combinemodels --------
#stick all models into a list
modellist = list(m1=m1eqs,m2a=m2aeqs,m3=m3eqs,m4=m4eqs)
#also make list for priors
priorlist = list(m1priors=m1priors,m2apriors=m2apriors,m3priors=m3priors,m4priors=m4priors)
# set up a list in which we'll store our results
fl = vector(mode = "list", length = length(modellist))



## ---- fittingsetup --------
#general settings for fitting
#you might want to adjust based on your computer
warmup = 6000
iter = warmup + floor(warmup/2)
max_td = 18 #tree depth
adapt_delta = 0.9999
chains = 5
cores  = chains
seed = 1234



## ---- startvalues --------
## Setting starting values
#starting values for model 1
startm1 = list(a0 = rep(2,Ntot), b0 = rep(0.5,Ntot), a1 = 0.5 , b1 = -0.5, sigma = 1)
#starting values for model 2a
startm2a = list(a0 = 2, b0 = 0.5, a1 = 0.5 , b1 = 0.5, sigma = 1)
#starting values for model 3
startm3 = startm1
#starting values for models 4
startm4 = list(mu_a = 2, sigma_a = 1, mu_b = 0, sigma_b = 1, a1 = 0.5 , b1 = -0.5, sigma = 1)
#put different starting values in list
#need to be in same order as models below
#one list for each chain, thus a 3-leveled list structure
#for each chain, we add jitter so they start at different values
startlist = list( rep(list(lapply(startm1,jitter,10)),chains),
                  rep(list(lapply(startm2a,jitter,10)),chains),
                  rep(list(lapply(startm3,jitter,10)),chains),
                  rep(list(lapply(startm4,jitter,10)),chains)
                  )


## ---- modelfitting --------
# fitting models
#loop over all models and fit them using ulam
for (n in 1:length(modellist))
{

  cat('************** \n')
  cat('starting model', names(modellist[n]), '\n')

  tstart=proc.time(); #capture current time

  fl[[n]]$fit <- brm(formula = modellist[[n]],
                   data = fitdat,
                   family = gaussian(),
                   prior = priorlist[[n]],
                   init = startlist[[n]],
                   control=list(adapt_delta=adapt_delta, max_treedepth = max_td),
                   sample_prior = TRUE,
                   chains=chains, cores = cores,
                   warmup = warmup, iter = iter,
                   seed = seed,
                   backend = "cmdstanr"
  )# end brm statement

  tend=proc.time(); #capture current time
  tdiff=tend-tstart;
  runtime_minutes=tdiff[[3]]/60;

  cat('model fit took this many minutes:', runtime_minutes, '\n')
  cat('************** \n')

  #add some more things to the fit object
  fl[[n]]$runtime = runtime_minutes
  fl[[n]]$model = names(modellist)[n]
}
# saving the results so we can use them later
filepath = fs::path("C:","Dropbox","datafiles","longitudinalbayes","brmsfits", ext="Rds")
#filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","brmsfits", ext="Rds")
saveRDS(fl,filepath)

## ---- additional-code -------
for (n in 1:length(fl)) {print(fl[[n]]$runtime)}

preprior1 <- get_prior(m1eqs,data=fitdat,family=gaussian())
preprior2 <- get_prior(m1eqs,data=fitdat,family=gaussian())




