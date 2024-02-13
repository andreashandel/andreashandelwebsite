# code for part 5 of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-5/

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

fitdat=list(id=simdat[[3]]$id,
            outcome = simdat[[3]]$outcome,
            time = simdat[[3]]$time,
            dose_adj = simdat[[3]]$dose_adj, 
            Nobs =  Nobs,
            Nind = Nind
            )



## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(here('posts','2024-02-15-longitudinal-multilevel-bayes-5',"stancode-2par.stan"), 
                                    pedantic=TRUE, 
                                    force_recompile=TRUE)


## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
#keeping values somewhat low to make things run reasonably fast
#for 'production' you would probably want to sample more 
#and set more stringent conditions
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
init_vals_1chain <- function() (list(mu_a = runif(1,1,3), 
                                     mu_b = runif(1,0,1),
                                     sigma_a = runif(1,1,10),
                                     sigma_b = runif(1,1,10),
                                     a1 = rnorm(1,-0.2,0.2),
                                     b1 = rnorm(1,-0.2,0.2),
                                     sigma = runif(1,1,10)))
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
print(res_m1$cmdstan_diagnose())


## ---- get_samples_m1 ----
#this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m1 <- res_m1$draws(inc_warmup = TRUE, format = "draws_df")


## ---- plot_par_m1 ----
# only main parameters, excluding parameters that we have for each individual, is too much
plotpars = c("a1","b1","a1_prior","b1_prior","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)

## ---- results_m1 ----
print(head(res_m1$summary(),15))
bp3 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
plot(bp3)
# just to get a picture that can be shown together with the post
#ggsave("featured.png",bp2)


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


## ---- make_predictions ----

#take originally saved data, which contains more variables 
# than the fitdat object, and process further
newdat <- simdat$m3
plotdat <- fitdat %>% data.frame()  %>%
           mutate(id = as.factor(id))  %>%
           mutate(dose = dose_cat)

#make new data for which we want predictions
#specifically, more time points so the curves are smoother
timevec = seq(from = 0.1, to = max(newdat$time), length=100)
Nind = max(newdat$id)
#new data used for predictions
preddat = data.frame( id = sort(rep(seq(1,Ntot),length(timevec))),
                      time = rep(timevec,Ntot),
                      dose_adj = 0
)
#add right dose information for each individual
for (k in 1:Nind)
{
  #dose for a given individual
  nowdose = unique(fitdat$dose_adj[fitdat$id == k])
  nowdose_cat = unique(fitdat$dose_cat[fitdat$id == k])
  #assign that dose
  #the categorical values are just for plotting
  preddat[(preddat$id == k),"dose_adj"] = nowdose
  preddat[(preddat$id == k),"dose_cat"] = nowdose_cat
}


# estimate and CI for parameter variation
# this is computed in the predicted-quantities block of the Stan code

post <- posterior::as_draws_array(fit2)

mu <- samp_m1 |>
  as_tibble() |>
  select(starts_with("virus")) |>
  apply(2, quantile, c(0.05, 0.5, 0.95)) |>
  t() %>%
  data.frame(time = newdat$time, .)  |> 
  tidyr::gather(pct, outcome, -time)


# estimate and CI for prediction intervals
# the predictions factor in additional uncertainty around the mean (mu)
# as indicated by sigma
# this is computed in the predicted-quantities block of the Stan code
outpred <- predict(nowmodel, newdata = preddat, probs = c(0.055, 0.945) )
  
# estimate and CI for parameter variation
# we ask for predictions for the new data generated above
linkmod <- rethinking::link(nowmodel, data = preddat)
  
  #computing mean and various credibility intervals
  #these choices are inspired by the Statistical Rethinking book
  #and purposefully do not include 95%
  #to minimize thoughts of statistical significance
  #significance is not applicable here since we are doing bayesian fitting
  modmean <- apply( linkmod$mu , 2 , mean )
  modPI79 <- apply( linkmod$mu , 2 , PI , prob=0.79 )
  modPI89 <- apply( linkmod$mu , 2 , PI , prob=0.89 )
  modPI97 <- apply( linkmod$mu , 2 , PI , prob=0.97 )
  
  # estimate and CI for prediction intervals
  # this uses the sim function from rethinking
  # the predictions factor in additional uncertainty around the mean (mu)
  # as indicated by sigma
  simmod <- rethinking::sim(nowmodel, data = preddat)
  
  # mean and credible intervals for outcome predictions
  # modmeansim should agree with above modmean values
  modmeansim <- apply( simmod , 2 , mean )
  modPIsim <- apply( simmod , 2 , PI , prob=0.89 )
  
  #place all predictions into a data frame
  fitpred = data.frame(id = as.factor(preddat$id),
                            dose = as.factor(preddat$dose_cat),
                            predtime = preddat$time,
                            Estimate = modmean,
                            Q79lo = modPI79[1,], Q79hi = modPI79[2,],
                            Q89lo = modPI89[1,], Q89hi = modPI89[2,],
                            Q97lo = modPI97[1,], Q97hi = modPI97[2,],
                            Qsimlo=modPIsim[1,], Qsimhi=modPIsim[2,]
  )



## ---- plot_predictions ----

predplot <- ggplot(data = fitpred, aes(x = predtime, y = Estimate, group = id, color = dose ) ) +
  geom_line() +
  #geom_ribbon( aes(x=time, ymin=Q79lo, ymax=Q79hi, fill = dose), alpha=0.6, show.legend = F) +
  geom_ribbon(aes(x=predtime, ymin=Q89lo, ymax=Q89hi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
  #geom_ribbon(aes(x=time, ymin=Q97lo, ymax=Q97hi, fill = dose), alpha=0.2, show.legend = F) +
  geom_ribbon(aes(x=predtime, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
  geom_point(data = plotdat, aes(x = time, y = outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Virus load",
       x = "days post infection") +
  theme_minimal() +
  ggtitle(title)
plot(predplot)

