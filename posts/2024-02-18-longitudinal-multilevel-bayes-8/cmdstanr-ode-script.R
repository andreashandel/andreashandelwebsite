# code for part 8 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-8/

## ---- packages --------
library('here') #for file loading
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting
library('bayesplot') #for plotting results


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
                 mu_e_mu = 1, mu_e_sd = 1,
                 a1_mu = 1, a1_sd = 0.2,
                 b1_mu = 1, b1_sd = 0.2,
                 g1_mu = 1, g1_sd = 0.2,
                 e1_mu = 1, e1_sd = 0.2,
                 V0_mu = 5, V0_sd = 1
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
            Ndose = length(unique(simdat[[3]]$dose_adj)),
            tstart = 0
            )
fitdat = c(fitdat,priorvals)

## ---- make_stanmodel -----
# make Stan model. 
stan_file <- here('posts','2024-02-18-longitudinal-multilevel-bayes-8',"stancode-ode.stan")
stanmod1 <- cmdstanr::cmdstan_model(stan_file = stan_file, 
                                    compile = TRUE,
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
             chains = 2,
             cores  = 2,
             seed = 1234,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(1234) #make inits reproducible
init_vals_1chain <- function() (list(mu_a = runif(1,2,4), 
                                     mu_b = runif(1,-1,-1),
                                     mu_g = runif(1,1,1),
                                     mu_e = runif(1,1,1),
                                     sigma_a = runif(1,0,1),
                                     sigma_b = runif(1,0,1),
                                     sigma_g = runif(1,0,1),
                                     sigma_e = runif(1,0,1),
                                     a0 = runif(Nind,3,1),
                                     b0 = runif(Nind,-1,1),
                                     g0 = runif(Nind,1,1),
                                     e0 = runif(Nind,0,1),
                                     a1 = runif(Ndose,0,1),
                                     b1 = runif(Ndose,0,1),
                                     g1 = runif(Ndose,0,1),
                                     e1 = runif(Ndose,0,1),
                                     sigma = runif(1,0,1),
                                     V0 = runif(Ndose,5,5)
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
                          init = fs_m1$init,
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
#filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode", ext="Rds")
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
# only picking a few parameters
plotpars = c("a0[1]","b0[1]","g0[1]","e0[1]","V0[1]","V0[2]","V0[3]","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
bp3 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)
plot(bp3)


## ---- prep_data_m1 ----
# data manipulation to get in shape for plotting
postdf <- samp_m1 %>% 
  select(!ends_with('prior')) %>% 
  select(!starts_with(".")) %>% 
  select(-"lp__") %>% 
  select(contains("[1]")) %>%
  rename_with(~ gsub("[1]", "", .x, fixed = TRUE) )
priordf <-  samp_m1 %>% 
  select(ends_with('prior')) %>% 
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE) ) 
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
# just to get a picture that can be shown together with the post
#ggsave("featured.png",m1_p1)



## ---- make_predictions ----
# averages and CI for 
# estimates of deterministic model trajectory
# for each observation
# this is computed in the transformed parameters block of the Stan code
mu <- samp_m1 |>
  select(starts_with("virus_pred")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() 
rownames(mu) <- NULL

# estimate and CI for prediction intervals
# the predictions factor in additional uncertainty around the mean (mu)
# as indicated by sigma
# this is computed in the predicted-quantities block of the Stan code
# the average of mu and preds should be more or less the same
# but preds will have wider uncertainty due to the residual variation sigma
preds <- samp_m1 |>
  select(starts_with("ypred")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() 
rownames(preds) <- NULL


## ---- plot_predictions ----

# change dose so it looks nicer in plot
dose = as.factor(fitdat$dose_adj)
levels(dose)[1] <- "low"
levels(dose)[2] <- "medium"
levels(dose)[3] <- "high"

#place everything into a data frame
fitpred = data.frame(id = as.factor(fitdat$id),
                     dose = dose,
                     time = fitdat$time,
                     Outcome = fitdat$outcome,
                     Estimate = mu[,2],
                     Qmulo = mu[,1], Qmuhi = mu[,3]
                     #Qsimlo = preds[,1], Qsimhi = preds[,3]
)

#make the plot
predplot <- ggplot(data = fitpred, aes(x = time, y = Estimate, group = id, color = dose ) ) +
  geom_line() +
  #geom_ribbon(aes(x=time, ymin=Qmulo, ymax=Qmuhi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
  #geom_ribbon(aes(x=time, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
  geom_point(aes(x = time, y = Outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Virus load",
       x = "days post infection") +
  theme_minimal() 
plot(predplot)






