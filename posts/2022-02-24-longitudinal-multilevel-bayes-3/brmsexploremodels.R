# exploration of fit objects code of brms tutorial (part 3)
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-3/


## ---- packages2 --------
library('dplyr') # for data manipulation
library('tidyr') # for data manipulation
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('brms') # for model fitting
library('posterior') #for post-processing
library('bayesplot') #for plots
library('fs') #for file path


## ---- loadfits --------
# loading list of previously saved fits.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local folder
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","brmsfits", ext="Rds")
if (!file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","brmsfits", ext="Rds")
}


fl <- readRDS(filepath)
# also load data file used for fitting
simdat <- readRDS("simdat.Rds")
#pull our the data set we used for fitting
#if you fit a different one of the simulated datasets, change accordingly
fitdat <- simdat$m3
#contains parameters used for fitting
pars <- simdat$m3pars


## ---- diagnostics ------
# Model 2a summary
#saving a bit of typing below
fit2 <- fl[[2]]$fit
summary(fit2)

## ---- traceplot ------
# Model 2a trace plots
plot(fit2)

## ---- traceplot-2 ------
# Another trace plot, using the bayesplot package
posterior <- rstan::extract(fit2$fit, inc_warmup = TRUE, permuted = FALSE)
bayesplot::mcmc_trace(posterior, n_warmup = 400, pars = variables(fit2)[c(1,2,3,4,5)])

## ---- trankplot ------
# Model 2a trank plots with bayesplot
bayesplot::mcmc_rank_overlay(fit2, pars = variables(fit2)[c(1,2,3,4,5)])

## ---- autocorrelationplot -----
bayesplot::mcmc_acf(fit2, pars = variables(fit2)[c(1,2,3,4,5)])


## ---- pairplot ------
# Model 2a pair plot
# Correlation between posterior samples of parameters
pairs(fit2)


## ---- mod_1_3_summary --------
#save some typing
fit1 <- fl[[1]]$fit
fit3 <- fl[[3]]$fit
summary(fit1)
summary(fit3)


## ---- mod_1_3_prior --------
#get priors and posteriors for models 1 and 3
m1prior <- prior_draws(fit1)
m1post <- as_draws_df(fit1)
m3prior <- prior_draws(fit3)
m3post <- as_draws_df(fit3)


## ---- mod_1_3_prior_plots --------
#showing density plots for a1

#make a data frame and get it in shape for ggplot
a1df <- data.frame(m1_prior = m1prior$b_alpha_dose_adj,
                   m1_post = m1post$b_alpha_dose_adj,
                   m3_prior = m3prior$b_alpha_dose_adj,
                   m3_post =  m3post$b_alpha_dose_adj) %>%
        pivot_longer(cols = everything(), names_to = c("model","type"), names_pattern = "(.*)_(.*)", values_to = "value")
# make plot
p1 <- a1df %>%
  ggplot() +
  geom_density(aes(x = value, color = model, linetype = type), size = 1) +
  theme_minimal()
plot(p1)
#save for display on post
ggsave(file = paste0("featured.png"), p1, dpi = 300, units = "in", width = 6, height = 6)


#showing density plots for b1
b1df <- data.frame(m1_prior = m1prior$b_beta_dose_adj,
                   m1_post = m1post$b_beta_dose_adj,
                   m3_prior = m3prior$b_beta_dose_adj,
                   m3_post =  m3post$b_beta_dose_adj) %>%
  pivot_longer(cols = everything(), names_to = c("model","type"), names_pattern = "(.*)_(.*)", values_to = "value")

p2 <- b1df %>%
  ggplot() +
  geom_density(aes(x = value, color = model, linetype = type), size = 1) +
  theme_minimal()
plot(p2)



## ---- mod_1_3_pair_plots --------
# a few parameters for each dose
#low dose
pairs(fit1, variable = variables(fit1)[c(1:4,25)])
#medium dose
pairs(fit1, variable = variables(fit1)[c(8:11,25)])
#high dose
pairs(fit1, variable = variables(fit1)[c(16:19,25)])





## ---- mod_1_3_pars --------
# model 1 first
fit1pars = posterior::summarize_draws(m1post, "mean", "sd", "quantile2", default_convergence_measures())

#only entries for the a0 parameters
a0post <- m1post %>% dplyr::select(starts_with('b_alpha_id'))
fit1a0mean <- mean(colMeans(a0post))
#only entries for the b0 parameters
b0post <- m1post %>% dplyr::select(starts_with('b_beta_id'))
fit1b0mean <- mean(colMeans(b0post))
fit1otherpars <- fit1pars %>% dplyr::filter(!grepl('_id',variable)) %>%
  dplyr::filter(!grepl('prior',variable))
print(fit1otherpars)
print(c(fit1a0mean,fit1b0mean))

# repeat for model 3
fit3pars = posterior::summarize_draws(m3post, "mean", "sd", "quantile2", default_convergence_measures())
#only entries for the a0 parameters
a0post <- m3post %>% dplyr::select(starts_with('b_alpha_id'))
fit3a0mean <- mean(colMeans(a0post))
#only entries for the b0 parameters
b0post <- m3post %>% dplyr::select(starts_with('b_beta_id'))
fit3b0mean <- mean(colMeans(b0post))
fit3otherpars <- fit3pars %>% dplyr::filter(!grepl('_id',variable)) %>%
  dplyr::filter(!grepl('prior',variable))
print(fit3otherpars)
print(c(fit3a0mean,fit3b0mean))


## ---- mod_1_3_comparison --------
fit13comp <- loo_compare(add_criterion(fit1,"waic"),
            add_criterion(fit3,"waic"),
            criterion = "waic")
print(fit13comp, simplify = FALSE)


## ---- mod_2a_exploration --------
m2post <- as_draws_df(fit2)
fit2pars = posterior::summarize_draws(m2post, "mean", "sd", "quantile2", default_convergence_measures())
fit2otherpars <- fit2pars %>% dplyr::filter(!grepl('prior',variable))
print(fit2otherpars)



## ---- mod_4_exploration --------
fit4 <- fl[[4]]$fit
m4prior <- prior_draws(fit4)
m4post <- as_draws_df(fit4)
summary(fit4)



## ---- mod_4_prior_plots --------
#showing density plots for a1 and b1
#make a data frame and get it in shape for ggplot
m4df <- data.frame(a1_prior = m4prior$b_alpha_dose_adj,
                   a1_post = m4post$b_alpha_dose_adj,
                   b1_prior = m4prior$b_beta_dose_adj,
                   b1_post = m4post$b_beta_dose_adj) %>%
  pivot_longer(cols = everything(), names_to = c("parameter","type"), names_pattern = "(.*)_(.*)", values_to = "value")
# make plot
p1 <- m4df %>%
  ggplot() +
  ylim(0, 10) + xlim(-2, 2) +
  geom_density(aes(x = value, color = parameter, linetype = type), adjust = 10, size = 1) +
  ggtitle('model 4, parameters a1 and b1') +
  theme_minimal()
plot(p1)


## ---- mod_4_pars --------
fit4pars = posterior::summarize_draws(m4post, "mean", "sd", "quantile2", default_convergence_measures())
fit4otherpars <- fit4pars %>% dplyr::filter(!grepl('_id',variable)) %>%
  dplyr::filter(!grepl('prior',variable)) %>%
  dplyr::filter(!grepl('z_',variable))

print(fit4otherpars)



## ---- mod_4_pair_plots --------
# a few parameters for each dose
#low dose
pairs(fit4, variable = variables(fit4)[c(1:4,25)])
#medium dose
pairs(fit4, variable = variables(fit4)[c(8:11,25)])
#high dose
pairs(fit4, variable = variables(fit4)[c(16:19,25)])




## ---- mod_all_comparison --------
fit1a <- add_criterion(fit1,c("waic","loo"))
fit2a <- add_criterion(fit2,c("waic","loo"))
fit3a <- add_criterion(fit3,c("waic","loo"))
fit4a <- add_criterion(fit4,c("waic","loo"))
compall1 <- loo_compare(fit1a,fit2a,fit3a,fit4a, criterion = "waic")
compall2 <- loo_compare(fit1a,fit2a,fit3a,fit4a, criterion = "loo")
print(compall1, simplify = FALSE)
print(compall2, simplify = FALSE)



## ---- priorexploration-1 --------
#defining model again
m2aeqs <- bf(outcome ~ exp(alpha)*log(time) - exp(beta)*time,
  alpha ~ 1 + dose_adj,
  beta  ~  1 + dose_adj,
  nl = TRUE)
preprior2 <- get_prior(m2aeqs,data=fitdat,family=gaussian())
postprior2 <- prior_summary(fit2)
print(preprior2)
print(postprior2)


## ---- priorexploration-2 --------
postprior1 <- prior_summary(fit1)
postprior3 <- prior_summary(fit3)
postprior4 <- prior_summary(fit4)
print(paste(nrow(postprior1),nrow(postprior3),nrow(postprior4)))


## ---- priorexploration-3 --------
names(m1post)

## ---- priorexploration-4 --------
print(postprior4)

## ---- priorexploration-5 --------
names(m4post)


## ---- computepredictions --------
#this will contain all the predictions from the different models
fitpred = vector(mode = "list", length = length(fl))

# load the data we used for fitting
simdat <- readRDS("simdat.Rds")
#pull our the data set we used for fitting
#if you fit a different one of the simulated datasets, change accordingly
fitdat <- simdat$m3
#small data adjustment for plotting
plotdat <- fitdat %>% data.frame() %>% mutate(id = as.factor(id)) %>% mutate(dose = dose_cat)


# we are looping over each fitted model
for (n in 1:length(fl))
{
  #get current model
  nowmodel = fl[[n]]$fit

  #make new data for which we want predictions
  #specifically, more time points so the curves are smoother
  timevec = seq(from = 0.1, to = max(fitdat$time), length=100)
  Ntot = max(fitdat$id)
  #data used for predictions
  preddat = data.frame( id = sort(rep(seq(1,Ntot),length(timevec))),
                        time = rep(timevec,Ntot),
                        dose_adj = 0
  )
  #add right dose information for each individual
  for (k in 1:Ntot)
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
  #brms equivalent to rethinking::link
  #doing 89% CI
  meanpred <- fitted(nowmodel, newdata = preddat, probs = c(0.055, 0.945) )

  # estimate and CI for prediction intervals
  # the predictions factor in additional uncertainty around the mean (mu)
  # as indicated by sigma
  # this is equivalent to rethinking::sim()
  outpred <- predict(nowmodel, newdata = preddat, probs = c(0.055, 0.945) )


  #place all predictions into a data frame
  #and store in a list for each model
  fitpred[[n]] = data.frame(id = as.factor(preddat$id),
                            dose = as.factor(preddat$dose_cat),
                            predtime = preddat$time,
                            Estimate = meanpred[,"Estimate"],
                            Q89lo = meanpred[,"Q5.5"],
                            Q89hi = meanpred[,"Q94.5"],
                            Qsimlo = outpred[,"Q5.5"],
                            Qsimhi = outpred[,"Q94.5"]
  )
}


#########################
# generate plots showing data and model predictions
#########################



## ---- makeplots --------
#storing all plots
plotlist = vector(mode = "list", length = length(fl))

#adding titles to plots
titles = c('model 1','model 2a','model 3','model 4')

#again looping over all models, making a plot for each
for (n in 1:length(fl))
{
  # ===============================================
  plotlist[[n]] <- ggplot(data = fitpred[[n]], aes(x = predtime, y = Estimate, group = id, color = dose ) ) +
    geom_line() +
    geom_ribbon(aes(x=predtime, ymin=Q89lo, ymax=Q89hi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
    geom_ribbon(aes(x=predtime, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
    geom_point(data = plotdat, aes(x = time, y = outcome, group = id, color = dose), shape = 1, size = 2) +
    scale_y_continuous(limits = c(-30,50)) +
    labs(y = "Virus load",
         x = "days post infection") +
    theme_minimal() +
    ggtitle(titles[n])
  ggsave(file = paste0(titles[n],".png"), plotlist[[n]], dpi = 300, units = "in", width = 7, height = 7)
}



#########################
# show the plots
#########################

## ---- showplots --------
plot(plotlist[[1]])
plot(plotlist[[3]])
plot(plotlist[[2]])
plot(plotlist[[4]])

