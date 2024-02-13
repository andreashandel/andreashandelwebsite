# code for part 6 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-6/

## ---- packages --------
library('here') #for file loading
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting
library('bayesplot') #for plotting results
library('loo') #for model diagnostics


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
# values for prior distributions
# allows for exploring different values without having to edit Stan model code
priorvals = list(mu_a_mu = 7, mu_a_sd = 3,
                 mu_b_mu = 1, mu_b_sd = 1,
                 mu_g_mu = 3, mu_g_sd = 1,
                 mu_e_mu = -3, mu_e_sd = 1
)

# all data as one list, this is how Stan needs it
fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            time = simdat[[3]]$time,
            dose_adj = simdat[[3]]$dose_adj[1:Nind], #first Nind values
            Nobs =  Nobs,
            Nind = Nind
            )
fitdat = c(fitdat,priorvals)

## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(here('posts','2024-02-16-longitudinal-multilevel-bayes-6',"stancode-4par.stan"), 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)


## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
fs_m1 = list(warmup = 1500,
             sampling = 1000, 
             max_td = 15, #tree depth
             adapt_delta = 0.999,
             chains = 5,
             cores  = 5,
             seed = 1234,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(1234) #make inits reproducible
init_vals_1chain <- function() (list(mu_a = runif(1,5,10), 
                                     mu_b = runif(1,1,2),
                                     mu_g = runif(1,1,4),
                                     mu_e = runif(1,-4,-2),
                                     sigma_a = runif(1,0,2),
                                     sigma_b = runif(1,0,2),
                                     sigma_g = runif(1,0,2),
                                     sigma_e = runif(1,0,2),
                                     a1 = rnorm(1,-0.1,0.1),
                                     b1 = rnorm(1,-0.1,0.1),
                                     g1 = rnorm(1,-0.1,0.1),
                                     e1 = rnorm(1,-0.1,0.1),
                                     sigma = runif(1,0,2)))
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

## ---- diagnose_m1 ----
res_m1$cmdstan_diagnose()

## ---- summarize_m1 ----
# uses posterior package 
print(head(res_m1$summary(),15))

## ---- get_samples_m1 ----
#this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m1 <- res_m1$draws(inc_warmup = TRUE, format = "draws_df")

## ---- plot_par_m1 ----
# only main parameters, excluding parameters that we have for each individual, is too much
plotpars = c("a1","b1","g1","e1","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
bp3 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)
plot(bp3)
# just to get a picture that can be shown together with the post
ggsave("featured.png",bp1)


## ---- prep_data_m1 ----
# data manipulation to get in shape for plotting
postdf <- samp_m1 %>% select(!ends_with('prior')) %>% select(!starts_with(".")) %>% select(-"lp__") %>% select(!contains("["))
priordf <- samp_m1 %>% select(ends_with('prior')) %>% rename_with(~ gsub("_prior", "", .x, fixed = TRUE) ) 
postlong <- tidyr::pivot_longer(data = postdf, cols = everything() , names_to = "parname", values_to = "value") %>% mutate(type = "posterior")
priorlong <- tidyr::pivot_longer(data = priordf, cols = everything() , names_to = "parname", values_to = "value") %>% mutate(type = "prior")
ppdf <- dplyr::bind_rows(postlong,priorlong)

## ---- prior_post_m1 ----
m1_p1 <- ppdf %>%
  ggplot() +
  geom_density(aes(x = value, color = type), linewidth = 1) +
  facet_wrap("parname", scales = "free") +
  theme_minimal()
plot(m1_p1)


## ---- obs_pred_m1 ----
ypred_df <- samp_m1 %>% select(starts_with("ypred"))
m1_p2 <- bayesplot::ppc_dens_overlay(fitdat$outcome, as.matrix(ypred_df))
plot(m1_p2)

## ---- loo_m1_part1 ----
# uses loo package 
loo_m1 <- res_m1$loo(cores = fs_m1$chains, save_psis = TRUE)
print(loo_m1)
plot(loo_m1)


## ---- loo_m1_part2 ----
ypred_df <- samp_m1 %>% select(starts_with("ypred"))
m1_p3 <- bayesplot::ppc_loo_pit_overlay(
  y = fitdat$outcome,
  yrep = as.matrix(ypred_df),
  lw = weights(loo_m1$psis_object)
)
plot(m1_p3)