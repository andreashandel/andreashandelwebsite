# code for part 6 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-6/

## ---- packages --------
library('here') #for file loading
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting
library('bayesplot') #for plotting results
library('posterior') #for post-processing
library('loo') #for model diagnostics


## ---- setup --------
############################################
# some general definitons and setup stuff
############################################
#setting random number seed for reproducibility
rngseed = 1234

# I'll be saving results so we can use them without always running the model
# Note that the files are often too large for standard Git/GitHub - where this project lives
# Git Large File Storage should be able to handle it
# I'm using a simple hack so I don't have to set up Git LFS
# I am saving these large file to a folder that is synced with Dropbox
# adjust accordingly for your setup
# filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr4par", ext="Rds")
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes")
filename = "cmdstanr4par.Rds"
stanfile <- here('posts','2024-02-16-longitudinal-multilevel-bayes-6',"stancode-4par.stan")


## ---- explore-model --------
# brief plotting of model to get idea for priors
t = seq(0.1,40,length=100) 
# all parameters are the log of their original values
alph = 30; # approximately the peak of virus
bet = 1.5; # approx. growth rate
gamm = 2; # approx. peak time
et = 0.3; # approx. decay rate
num  = 2*exp(alph)
d1 = exp( - exp(bet)*(t - exp(gamm)) )
d2 =  exp( exp(et) * (t - exp(gamm)) )
mu = log( num /(d1 + d2) ) 
plot(t,mu, type = "l") #looks somewhat like virus load in acute 


## ---- data --------
# adjust as necessary
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
Nind = length(unique(simdat$m3$id))
Ntot =  length(simdat$m3$id)
# values for prior distributions
# allows for exploring different values without having to edit Stan model code
priorvals = list(mu_a_mu = 30, mu_a_sd = 5,
                 mu_b_mu = 1.5, mu_b_sd = 1,
                 mu_g_mu = 2, mu_g_sd = 1,
                 mu_e_mu = 0.5, mu_e_sd = 1,
                 a1_mu = 0.5, a1_sd = 1,
                 b1_mu = 0.1, b1_sd = 1,
                 g1_mu = 0.1, g1_sd = 1,
                 e1_mu = -0.1, e1_sd = 1
)

# all data as one list, this is how Stan needs it
fitdatbase=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            time = simdat[[3]]$time,
            dose_adj = simdat[[3]]$dose_adj[1:Nind], #first Nind values
            Ntot =  Ntot,
            Nind = Nind
            )
fitdat = c(fitdatbase,priorvals)

## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(stanfile, 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)


## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
fs_m1 = list(warmup = 1500,
             sampling = 1000, 
             max_td = 20, #tree depth
             adapt_delta = 0.99999,
             chains = 5,
             cores  = 5,
             seed = rngseed,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(rngseed) #make inits reproducible
init_vals_1chain <- function() (list(mu_a = runif(1,25,35), 
                                     mu_b = runif(1,1,2),
                                     mu_g = runif(1,1.5,2.5),
                                     mu_e = runif(1,0,1),
                                     sigma_a = runif(1,0,1),
                                     sigma_b = runif(1,0,1),
                                     sigma_g = runif(1,0,1),
                                     sigma_e = runif(1,0,1),
                                     a0 = runif(Nind,25,35),
                                     b0 = runif(Nind,1,2),
                                     g0 = runif(Nind,1.5,2.5),
                                     e0 = runif(Nind,0,1),
                                     a1 = runif(1,0.5,0.6),
                                     b1 = runif(1,0.1,0.1),
                                     g1 = runif(1,0.1,0.1),
                                     e1 = runif(1,-0.1,-0.1),
                                     sigma = runif(1,0,1))
                                )
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
                          adapt_delta = fs_m1$adapt_delta,
                          output_dir = filepath
)

## ---- savefits ----
res_m1$save_object(fs::path(filepath,filename))


## ---- loadfits --------
# loading previously saved fit.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local cloud-synced folder
# adjust accordingly for your setup
res_m1 <- readRDS(fs::path(filepath,filename))


## ---- diagnose_m1 ----
res_m1$cmdstan_diagnose()

## ---- get_samples_m1 ----
#this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m1 <- res_m1$draws(inc_warmup = TRUE, format = "draws_df")

## ---- plot_par_m1 ----
# only main parameters, excluding parameters that we have for each individual, is too much
plotpars = c("a1","b1","g1","e1","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)


## ---- results_m1 ----
print(head(res_m1$summary(),15))
bp3 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
plot(bp3)


## ---- prep_data_m1 ----
# data manipulation to get in shape for plotting
postdf1 <- samp_m1 %>% 
  select(!ends_with('prior')) %>% 
  select(!starts_with(".")) %>% 
  select(-"lp__") %>% 
  select(!contains("[")) 
# awkward way of getting some further parameters
# namely values from first individual for a0,b0,g0,e0
postdf2 <- samp_m1 %>%
           select(contains("0[1]")) %>%
           rename_with(~ gsub("[1]", "", .x, fixed = TRUE) )
postdf <- cbind(postdf1, postdf2) 
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


## ---- obs_pred_m1 ----
ypred_df <- samp_m1 %>% select(starts_with("ypred"))
m1_p2 <- bayesplot::ppc_dens_overlay(fitdat$outcome, as.matrix(ypred_df))
plot(m1_p2)
# just to get a picture that can be shown together with the post
#ggsave("featured.png",m1_p2)


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
                     Qmulo = mu[,1], Qmuhi = mu[,3],
                     Qsimlo = preds[,1], Qsimhi = preds[,3]
)

#make the plot
predplot <- ggplot(data = fitpred, aes(x = time, y = Estimate, group = id, color = dose ) ) +
  geom_line() +
  geom_ribbon(aes(x=time, ymin=Qmulo, ymax=Qmuhi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
  #geom_ribbon(aes(x=time, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
  geom_point(aes(x = time, y = Outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Virus load",
       x = "days post infection") +
  theme_minimal() 
plot(predplot)



## ---- data_m2 --------
# need to update priors
priorvals2 = list(a0_mu = 30, a0_sd = 5,
                 b0_mu = 1.5, b0_sd = 1,
                 g0_mu = 2, g0_sd = 1,
                 e0_mu = 0.5, e0_sd = 1
                )
fitdat2 = c(fitdatbase,priorvals2)

## ---- initialconditions_m2 ----
# need different initial values
init_vals_1chain <- function() (list(a0_mu = runif(Nind,25,35), 
                                    b0_mu = runif(Nind,1,2),
                                     g0_mu = runif(Nind,1.5,2.5),
                                     e0_mu = runif(Nind,0,1),
                                   
                                     a0_sd = runif(Nind,5,5),
                                     b0_sd = runif(Nind,1,1),
                                     g0_sd = runif(Nind,1,1),
                                     e0_sd = runif(Nind,1,1),
                                     sigma = runif(1,0,1))
                                )
inits = NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m1$init = inits

## ---- make_stanmodel2 -----
# make Stan model. 
stanfile <- here('posts','2024-02-16-longitudinal-multilevel-bayes-6',"stancode-4par-simple.stan")
stanmod2 <- cmdstanr::cmdstan_model(stanfile, 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)




## ---- run_m2 ----
res_m2 <- stanmod2$sample(data = fitdat2,
                          chains = fs_m1$chains,
                          init = fs_m1$init,
                          seed = fs_m1$seed,
                          parallel_chains = fs_m1$chains,
                          iter_warmup = fs_m1$warmup,
                          iter_sampling = fs_m1$sampling,
                          save_warmup = fs_m1$save_warmup,
                          max_treedepth = fs_m1$max_td,
                          adapt_delta = fs_m1$adapt_delta,
                          output_dir = filepath
)

## ---- savefits2 ----
filename = "cmdstanr4par-simple.Rds"
res_m2$save_object(fs::path(filepath,filename))

## ---- loadfits2 --------
res_m2 <- readRDS(fs::path(filepath,filename))

## ---- get_samples_m2 ----
#this uses the posterior package to get draws
samp_m2 <- res_m2$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m2 <- res_m2$draws(inc_warmup = TRUE, format = "draws_df")



## ---- make_predictions_m2 ----
mu2 <- samp_m2 |>
  select(starts_with("virus_pred")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() 
rownames(mu) <- NULL
preds2 <- samp_m2 |>
  select(starts_with("ypred")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() 
rownames(preds) <- NULL


## ---- plot_predictions_m2 ----
# change dose so it looks nicer in plot
dose = as.factor(fitdat2$dose_adj)
levels(dose)[1] <- "low"
levels(dose)[2] <- "medium"
levels(dose)[3] <- "high"
fitpred2 = data.frame(id = as.factor(fitdat2$id),
                     dose = dose,
                     time = fitdat2$time,
                     Outcome = fitdat2$outcome,
                     Estimate = mu[,2],
                     Qmulo = mu2[,1], Qmuhi = mu2[,3],
                     Qsimlo = preds2[,1], Qsimhi = preds2[,3]
)

#make the plot
predplot2 <- ggplot(data = fitpred2, aes(x = time, y = Estimate, group = id, color = dose ) ) +
  geom_line() +
  #geom_ribbon(aes(x=time, ymin=Qmulo, ymax=Qmuhi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
  #geom_ribbon(aes(x=time, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
  geom_point(aes(x = time, y = Outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Virus load",
       x = "days post infection") +
  theme_minimal() 
plot(predplot2)

