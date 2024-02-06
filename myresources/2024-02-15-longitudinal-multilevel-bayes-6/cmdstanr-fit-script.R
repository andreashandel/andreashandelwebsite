# code for part 5 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-5/

## ---- packages --------
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting

## ---- data --------
# adjust as necessary
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
# using dataset 3 for fitting
# also removing anything in the dataframe that's not used for fitting
# makes the fitting more robust
# We format the data slightly differently compared to the prior examples
# We only store the dose for each individual
# And also add number of observations and individuals
Nind = length(unique(simdat$m3$id))
Nobs =  length(simdat$m3$id)

fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            time = simdat[[3]]$time,
            dose_adj = simdat[[3]]$dose_adj[1:Nind], #first Nind values
            Nobs =  Nobs,
            Nind = Nind
            )



## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model("stancode-2par.stan", 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)


## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
fs_m1 = list(warmup = 1000,
             sampling = 1000, 
             max_td = 15, #tree depth
             adapt_delta = 0.999,
             chains = 4,
             cores  = 4,
             seed = 1234,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(1234) #make inits reproducible
init_vals_1chain <- function() (list(m = runif(1,90,150), 
                                     k = runif(1,1,2), 
                                     sigma = runif(1,8,12)))
inits = NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m1$init = inits


## ---- run_m1 ----
res_m1 <- stanmod1$sample(data = fitdat,
                          chains = fs_m1$chains,
                          init = fs_m1$init,
                          seed = fs_m1$seed,
                          parallel_chains = fs_m1$chains,
                          iter_warmup = fs_m1$warmup,
                          iter_sampling = fs_m1$sampling,
                          save_warmup = fs_m1$save_warmup,
                          max_treedepth = fs_m1$max_td,
                          adapt_delta = fs_m1$adapt_delta
)





# saving the list of results so we can use them later
# the file is too large for standard Git/GitHub
# Git Large File Storage should be able to handle it
# I'm using a simple hack so I don't have to set up Git LFS
# I am saving these large file to a folder that is synced with Dropbox
# adjust accordingly for your setup
#filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
# filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","cmdstanrfits", ext="Rds")
# saveRDS(fl,filepath)


