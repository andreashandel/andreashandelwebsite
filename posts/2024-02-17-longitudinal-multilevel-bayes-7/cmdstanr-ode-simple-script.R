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
    dU = - bet*y[1]*y[3]
    dI = bet*y[1]*y[3] - gamm*y[2]
    dV = alph*y[2] - et*y[3]
    return(list(c(dU,dI,dV)))
}
# initial conditions, on log scale 
# Uninfected cells, infected cells, virus
# note that I'm starting with 1 infected cells so I can take the log
# also, the starting value for V has the extra exponent because I'll be
# fitting this as a free parameter and therefore also need to do to 
# exponentiation trick
y0 = log(c(1e8,1,exp(5)))
t = seq(0,10,length=100)
# some parameter values that lead to a decent curve
alph = 3
bet = -1
gamm = 1 
et = 1
# all parameters
parms = exp(c(alph=alph,bet=bet,gamm=gamm,et=et))
# reproductive number, gives an indication of virus growth
# needs to be >1 
R0 = exp(bet)*exp(alph)*y0[3]/(exp(gamm)*exp(et))
print(R0)
# run the ODE model
oderes <- ode(y = y0, t = t, func = odemod, parms = parms)
# look at virus load trajectory
# note that it's already on the log scale
plot(oderes[,1],oderes[,4])


## ---- data --------
# adjust as necessary
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
# using dataset 3 for fitting
# formatted for Stan use
Ntot =  length(simdat$m3$id) #total observations
Nind = length(unique(simdat$m3$id)) #number of individuals
Nobs = as.numeric(table(simdat[[3]]$id)) #number of observations per individual, the same for everyone here
Ndose = length(unique(simdat$m3$dose_adj))
# values for prior distributions
# allows for exploring different values without having to edit Stan model code
priorvals = list(a0_mu = 3, a0_sd = 1,
                 b0_mu = -1, b0_sd = 1,
                 g0_mu = 1, g0_sd = 1,
                 e0_mu = 1, e0_sd = 1,
                 V0_mu = exp(5), V0_sd = 10
)

# all data as one list, this is how Stan needs it
fitdat=list(id=simdat[[3]]$id, #an ID for each individual, for indexing
            outcome = simdat[[3]]$outcome, #all outcomes
            time = simdat[[3]]$time, #all times
            dose_level = as.numeric(factor(simdat[[3]]$dose_adj[1:Nind])), #dose category for each individual
            Ntot =  Ntot,
            Nobs =  Nobs,
            Nind = Nind,
            Ndose = Ndose,
            tstart = 0
            )
fitdat = c(fitdat,priorvals)


## ---- make_stanmodel -----
# make Stan model. 
stan_file <- here('posts','2024-02-17-longitudinal-multilevel-bayes-7',"stancode-ode-simple.stan")
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
             chains = 1,
             cores  = 1,
             seed = 1234,
             save_warmup = TRUE)


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# a different sample will be drawn for each chain
# there's probably a better way to do that than a for loop
set.seed(1234) #make inits reproducible
init_vals_1chain <- function() (list(a0 = runif(Nind,3,4),
                                     b0 = runif(Nind,-1,-1),
                                     g0 = runif(Nind,1,1),
                                     e0 = runif(Nind,1,1),
                                     sigma = runif(1,0,1),
                                     V0 = runif(Ndose,4,6)
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
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode-simple", ext="Rds")
#filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode-simple", ext="Rds")
res_m1$save_object(file=filepath)


## ---- loadfits --------
# loading previously saved fit.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local cloud-synced folder
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode-simple", ext="Rds")
if (!fs::file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","cmdstanr-ode-simple", ext="Rds")
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
# only main parameters, excluding parameters that we have for each individual, is too much
plotpars = c("a1","b1","g1","e1","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp2 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
bp3 <- bayesplot::mcmc_pairs(samp_m1, pars = plotpars)
plot(bp1)
plot(bp2)
plot(bp3)


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




