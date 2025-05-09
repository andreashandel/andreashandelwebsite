# code for part 5 of longitudinal Bayesian fitting tutorial
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-5/

## ---- packages --------
library("here") # for file loading
library("dplyr") # for data manipulation
library("ggplot2") # for plotting
library("fs") # for file path
library("cmdstanr") # for model fitting
library("bayesplot") # for plotting results
library("posterior") # for post-processing
library("loo") # for model diagnostics


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
filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes")
#filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes")
filename = "cmdstanr2par.Rds"
stanfile <- here('posts','2024-02-15-longitudinal-multilevel-bayes-5',"stancode-2par.stan")



## ---- data --------
# adjust as necessary
simdatloc <- here::here("posts", "2022-02-22-longitudinal-multilevel-bayes-1", "simdat.Rds")
simdat <- readRDS(simdatloc)
# again using dataset 3 for fitting
# some formatting to get it into the shape needed for cmdstanr
Nind <- length(unique(simdat$m3$id)) # number of individuals
Nobs <- length(simdat$m3$id) # total number of observations
fitdat <- list(
  id = simdat[[3]]$id,
  outcome = simdat[[3]]$outcome,
  time = simdat[[3]]$time,
  dose_adj = simdat[[3]]$dose_adj,
  Nobs = Nobs,
  Nind = Nind
)

## ---- show_stancode -----
# not done this way currently, using alternative, see Quarto file
#print(stanmod1)

## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(stanfile, 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)



## ---- fitconditions ----
# settings for fitting
# keeping values somewhat low to make things run reasonably fast
# for 'production' you would probably want to sample more
# and set more stringent conditions
fs_m1 <- list(
  warmup = 1500,
  sampling = 2000,
  max_td = 18, # tree depth
  adapt_delta = 0.9999,
  chains = 5,
  cores = 5,
  seed = rngseed,
  save_warmup = TRUE
)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(rngseed) #make inits reproducible
init_vals_1chain <- function() {
  (list(
    mu_a = runif(1, 1, 3),
    mu_b = runif(1, 0, 1),
    sigma_a = runif(1, 1, 10),
    sigma_b = runif(1, 1, 10),
    a1 = rnorm(1, -0.2, 0.2),
    b1 = rnorm(1, -0.2, 0.2),
    sigma = runif(1, 1, 10)
  ))
}
inits <- NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] <- init_vals_1chain()
}
fs_m1$init <- inits


## ---- run_m1 ----
res_m1 <- stanmod1$sample(
  data = fitdat,
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
print(res_m1$cmdstan_diagnose())


## ---- get_samples_m1 ----
# this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m1 <- res_m1$draws(inc_warmup = TRUE, format = "draws_df")


## ---- plot_par_m1 ----
# only main parameters 
# excluding parameters that we have for each individual, is too much
plotpars <- c("a1", "b1", "a1_prior", "b1_prior", "sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)
# just to get a picture that can be shown together with the post
# ggsave("featured.png",bp2)

## ---- results_m1 ----
print(head(res_m1$summary(), 15))
bp3 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
plot(bp3)


## ---- prep_data_m1 ----
# data manipulation to get in shape for plotting
postdf <- samp_m1 %>%
  select(!ends_with("prior")) %>%
  select(!starts_with(".")) %>%
  select(-"lp__") %>%
  select(!contains("["))
priordf <- samp_m1 %>%
  select(ends_with("prior")) %>%
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE))
postlong <- tidyr::pivot_longer(data = postdf, cols = everything(), names_to = "parname", values_to = "value") %>% mutate(type = "posterior")
priorlong <- tidyr::pivot_longer(data = priordf, cols = everything(), names_to = "parname", values_to = "value") %>% mutate(type = "prior")
ppdf <- dplyr::bind_rows(postlong, priorlong)

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
dose <- as.factor(fitdat$dose_adj)
levels(dose)[1] <- "low"
levels(dose)[2] <- "medium"
levels(dose)[3] <- "high"

# place everything into a data frame
fitpred <- data.frame(
  id = as.factor(fitdat$id),
  dose = dose,
  time = fitdat$time,
  Outcome = fitdat$outcome,
  Estimate = mu[, 2],
  Qmulo = mu[, 1], Qmuhi = mu[, 3],
  Qsimlo = preds[, 1], Qsimhi = preds[, 3]
)

# make the plot
predplot <- ggplot(data = fitpred, aes(x = time, y = Estimate, group = id, color = dose)) +
  geom_line() +
  geom_ribbon(aes(x = time, ymin = Qmulo, ymax = Qmuhi, fill = dose, color = NULL), alpha = 0.3, show.legend = F) +
  geom_ribbon(aes(x = time, ymin = Qsimlo, ymax = Qsimhi, fill = dose, color = NULL), alpha = 0.1, show.legend = F) +
  geom_point(aes(x = time, y = Outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
  scale_y_continuous(limits = c(-30, 50)) +
  labs(
    y = "Virus load",
    x = "days post infection"
  ) +
  theme_minimal()
plot(predplot)
