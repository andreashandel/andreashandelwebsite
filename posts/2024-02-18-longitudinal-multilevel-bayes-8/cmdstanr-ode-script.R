# code for part 7 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-7/

## ---- packages --------
library('here') #for file loading
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting
library('bayesplot') #for plotting results
library('loo') #for model diagnostics
library("deSolve") #to explore the ODE model in R

## ---- explore-model --------
# brief plotting of model to get idea for priors
odemod <- function(t,y,parms) 
  {
    alph = parms[1]; bet = parms[2]; gamm = parms[3]; et = parms[4]
    dU = -exp(bet)*y[1]*y[3]
    dI = exp(bet)*y[1]*y[3] - exp(gamm)*y[2]
    dV = exp(alph)*y[2] - exp(et) * y[3]
    return(list(c(dU,dI,dV)))
}
y0 = log(c(1e8,1,100))
t = seq(0,10,length=100)
alph = 3
bet = -1
gamm = 1 
et = 1
parms = c(alph=alph,bet=bet,gamm=gamm,et=et)
R0 = exp(bet)*exp(alph)*y0[3]/(exp(gamm)*exp(et))
print(R0)
oderes <- ode(y = y0, t = t, func = odemod, parms = parms)
plot(oderes[,1],(oderes[,4]))


## ---- data --------
# adjust as necessary
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
# using dataset 3 for fitting
# formatted for Stan use
Ntot =  length(simdat$m3$id) #total observations
Nind = length(unique(simdat$m3$id)) #number of individuals
Nobs = as.numeric(table(simdat[[3]]$id)) #number of observations per individual, the same for everyone here
# values for prior distributions
# allows for exploring different values without having to edit Stan model code
priorvals = list(mu_a_mu = 10, mu_a_sd = 1,
                 mu_b_mu = -23, mu_b_sd = 1,
                 mu_g_mu = 1, mu_g_sd = 1,
                 mu_e_mu = 1, mu_e_sd = 1
)

# all data as one list, this is how Stan needs it
fitdat=list(id=simdat[[3]]$id, #an ID for each individual, for indexing
            outcome = simdat[[3]]$outcome, #all outcomes
            time = simdat[[3]]$time, #all times
            dose_adj = simdat[[3]]$dose_adj[1:Nind], #first Nind values
            dose_level = as.numeric(factor(simdat[[3]]$dose_adj[1:Nind])), #dose category for each individual
            Ntot =  Ntot,
            Nobs =  Nobs,
            Nind = Nind,
            tstart = 0
            )
fitdat = c(fitdat,priorvals)

## ---- make_stanmodel-simple -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(here('posts','2024-02-17-longitudinal-multilevel-bayes-7',"stancode-ode-simple.stan"), 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)



## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(here('posts','2024-02-17-longitudinal-multilevel-bayes-7',"stancode-ode.stan"), 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)


## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
fs_m1 = list(warmup = 1500,
             sampling = 1000, 
             max_td = 15, #tree depth
             adapt_delta = 0.99,
             chains = 1,
             cores  = 1,
             seed = 1234,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(1234) #make inits reproducible
init_vals_1chain <- function() (list(mu_a = runif(1,9,11), 
                                     mu_b = runif(1,-24,-22),
                                     mu_g = runif(1,1,2),
                                     mu_e = runif(1,1,2),
                                     sigma_a = runif(1,0,1),
                                     sigma_b = runif(1,0,1),
                                     sigma_g = runif(1,0,1),
                                     sigma_e = runif(1,0,1),
                                     a0 = runif(Nind,-1,1),
                                     b0 = runif(Nind,-1,1),
                                     g0 = runif(Nind,-1,1),
                                     e0 = runif(Nind,-1,1),
                                     a1 = runif(1,3,3),
                                     b1 = runif(1,-1,-1),
                                     g1 = runif(1,1,1),
                                     e1 = runif(1,1,1),
                                     sigma = runif(1,0,1),
                                     V0 = runif(3,5,5)
                                     ))

inits = NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m1$init = inits


## ---- run_m1 ----
res_m1 <- stanmod1$sample(data = fitdat,
                          chains = fs_m1$chains,
                          #init = fs_m1$init,
                          seed = fs_m1$seed,
                          parallel_chains = fs_m1$chains,
                          iter_warmup = fs_m1$warmup,
                          iter_sampling = fs_m1$sampling,
                          save_warmup = fs_m1$save_warmup,
                          max_treedepth = fs_m1$max_td,
                          adapt_delta = fs_m1$adapt_delta
)

## ---- savefits ----
# saving the list of results so we can use them later
# the file is too large for standard Git/GitHub
# Git Large File Storage should be able to handle it
# I'm using a simple hack so I don't have to set up Git LFS
# I am saving these large file to a folder that is synced with Dropbox
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode", ext="Rds")
if (!fs::file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode", ext="Rds")
}
res_m1$save_object(file=filepath)


## ---- loadfits --------
# loading previously saved fit.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local cloud-synced folder
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode", ext="Rds")
if (!fs::file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode", ext="Rds")
}
res_m1 <- readRDS(filepath)

