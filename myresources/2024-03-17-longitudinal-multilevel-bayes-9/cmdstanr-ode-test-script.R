# code for part XX of longitudinal Bayesian fitting tutorial 
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-XX/

## ---- packages --------
library('here') #for file loading
library('dplyr') # for data manipulation
library('tidyr') # for data manipulation
library('ggplot2') # for plotting
library('fs') #for file path
library('cmdstanr') #for model fitting
library("deSolve") #to explore the ODE model in R

## ---- setup --------
############################################
# some general definitons and setup stuff
############################################
#setting random number seed for reproducibility
rngseed = 1234
# File locations and names
# adjust as needed
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes")
#filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes")
filename = "cmdstanr-ode-testing.Rds"
stanfile1 <- here('myresources','2024-03-17-longitudinal-multilevel-bayes-9',"stancode-ode-testing1.stan")
stanfile2 <- here('myresources','2024-03-17-longitudinal-multilevel-bayes-9',"stancode-ode-testing2.stan")
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')

## ---- data --------
# loading data, path set in setup
simdat <- readRDS(simdatloc)
# using dataset 3 for fitting
# need to reformat all data as one list, this is how Stan needs it
dat <- simdat[[3]] %>% dplyr::arrange(id,time)

Ntot =  length(dat$id) #total observations
Nind = length(unique(dat$id)) #number of individuals
Nobs = as.numeric(table(dat$id)) #number of observations per individual, the same for everyone here
Ndose = length(unique(dat$dose_adj))
# since Stan needs the data as one long vector
# we make two vectors that indicate for each individual where the first
# and last observation is
start = rep(0,Nind)
stop = start
start[1] = 1
stop[1] = Nobs[1]
for (i in 2:Nind) {
  start[i] = start[i - 1] + Nobs[i - 1];
  stop[i] = stop[i - 1] + Nobs[i];
}
fitdat=list(id=dat$id, #an ID for each individual, for indexing
            outcome = dat$outcome, #all outcomes
            time = dat$time, #all times
            dose_adj = dat$dose_adj,
            dose_level = as.numeric(factor(dat$dose_adj)), 
            Ntot =  Ntot,
            Nobs =  Nobs,
            Nind = Nind,
            Ndose = Ndose,
            start = start,
            stop = stop,
            tstart = 0
)


## ---- odemodel1 --------
# regular model
odemod1 <- function(t,y,parms) 
  {
    alph = parms[1]; bet = parms[2]; gamm = parms[3]; et = parms[4]

    dU = - bet*y[1]*y[3]
    dI = bet*y[1]*y[3] - gamm*y[2]
    dV = alph*y[2] - et*y[3]
    return(list(c(dU,dI,dV)))
}

## ---- odemodel2 --------
# simulating in log space
odemod2 <- function(t,y,parms) 
{
  alph = parms[1]; bet = parms[2]; gamm = parms[3]; et = parms[4]
  
  dlogU = - bet*exp(y[3])
  dlogI = bet*exp(y[1] + y[3] - y[2]) - gamm
  dlogV = alph*exp(y[2] - y[3]) - et
  return(list(c(dlogU,dlogI,dlogV)))
}


## ---- priors --------
# values for prior distributions
pv = list(a0_mu = 1, a0_sd = 0.1,
                 b0_mu = -1, b0_sd = 0.1,
                 g0_mu = 0, g0_sd = 1,
                 e0_mu = 1, e0_sd = 1
)


## ---- r-simulation --------
Nsamp = 20 #number of samples for which to generate trajectories
timevec = seq(0,max(dat$time),length=1000) #time
odeall1 = NULL #set up empty array to hold all simulations
odeall2 = NULL #same empty array for model version 2
U0 = 1e12; #number of uninfected cells
I0 = 1e-15; #set to small non-zero so we can take log
V0 = 1e-3; #same value for everyone - doesn't agree with data but easy to try
y0 = c(U0,I0,V0)  
# scaling such that parameter distributions are close to 1
as = 12; #scaling of a0
bs = 35; #scaling of b0
set.seed(rngseed) #set seed for reproducibility

# run both models for the indicated number of samples
# could be done more efficiently without a loop
for (n in 1:Nsamp)
{
  # draw samples for parameters
  alph = exp(as*rnorm(1,pv$a0_mu, pv$a0_sd))
  bet = exp(bs*rnorm(1,pv$b0_mu, pv$b0_sd))
  gamm = exp(rnorm(1,pv$g0_mu, pv$g0_sd)) 
  et = exp(rnorm(1,pv$e0_mu, pv$e0_sd))
  odeparms = c(alph=alph,bet=bet,gamm=gamm,et=et)  
  # run the ODE model twice, in linear and log space
  # should give the same results 
  oderes1 <- ode(y = y0, t = timevec, func = odemod1, parms = odeparms, method = "bdf")
  oderes2 <- ode(y = log(y0), t = timevec, func = odemod2, parms = odeparms, method = "bdf")
  # combine all results, also add individual ID as n
  odeall1=rbind(odeall1, cbind(oderes1,n))
  odeall2=rbind(odeall2, cbind(oderes2,n))
}

## ---- r-simulation-plots --------
#convert to data frame
odeall1 = data.frame(odeall1)
# make plot showing all trajectories and the data
p1 <- odeall1 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id)), data = dat) +
  #scale_y_log10() + 
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p1)
# repeat for model 2
# result should look same as for model 1
odeall2 = data.frame(odeall2)
p2 <- odeall2 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=X3, group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p2)





## ---- make-stanmodel1 -----
# make Stan model 
stanmod1 <- cmdstanr::cmdstan_model(stan_file = stanfile1, 
                                    compile = TRUE,
                                    pedantic=TRUE,
                                    force_recompile=TRUE)


## ---- sim-conditions ----
#settings specific for simulation, no fitting is done here 
fs_m1 = list(warmup = 1,
             sampling = 40, 
             chains = 1,
             cores  = 1,
             seed = 1234
             )


## ---- initialconditions ----
# separate definition of initial values, added to fs_m1 structure 
# fixing here to the mean of the distributions
init_vals_1chain <- function() (list(a0 = rep(pv$a0_mu,Nind),
                                     b0 = rep(pv$b0_mu,Nind),
                                     g0 = rep(pv$g0_mu,Nind),
                                     e0 = rep(pv$e0_mu,Nind),
                                     sigma = 1
))
inits = NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m1$init = inits


## ---- run-code1-m1 ----
# adding priors to data
fitdatstan = c(fitdat,
               pv,
               as = 12, #scaling of a0
               bs = 35 #scaling of b0
              )

# run the model to generate simulations
res_m1 <- stanmod1$sample(data = fitdatstan,
                          chains = fs_m1$chains,
                          init = fs_m1$init,
                          seed = fs_m1$seed,
                          parallel_chains = fs_m1$chains,
                          iter_warmup = fs_m1$warmup,
                          iter_sampling = fs_m1$sampling,
                          output_dir = filepath,
                          fixed_param = TRUE
                          )




## ---- get_samples_m1 ----
#this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")


## ---- par_test_stan ----
#get samples of parameters for first individual 
x1 <- samp_m1 %>% select(c("alph[1]","bet[1]","gamm[1]","et[1]")) 
odeall3 = NULL
for (n in 1:nrow(x1)){
  # virus load varies by dose, but barely so keeping it fixed here
  parms = c(alph=x1$`alph[1]`[n],bet=x1$`bet[1]`[n],gamm=x1$`gamm[1]`[n],et=x1$`et[1]`[n])
  # run model for each set of parameters for each individual
  oderes3 <- ode(y = y0, t = timevec, func = odemod1, parms = parms)
  # combine all results, also add individual ID as n and add dose for a given individual
  odeall3 = rbind(odeall3, cbind(oderes3, n, dat$dose_cat[match(n,dat$id)]))
}
#plot virus load curve for each individual from ODE model
odeall3 = data.frame(odeall3)


## ---- pred_test_stan ----
#get predictions for first model  
x2 <- samp_m1 %>% select( contains("virus_pred1")) 
#observations of all samples for first individual only
x3 <- data.frame(x2[,1:Nobs[1]]) 
# turn into a long data frame that contains value for each individual and time
# some transformation shennanigans to get the values in the right order
df4 <- data.frame(time = rep(dat$time, nrow(samp_m1)), 
                 sample = rep(1:nrow(samp_m1), each = Nobs[1]),
                 value = as.vector(t(as.matrix(x3)))
                 )


## ---- fig_test_stan ----
# plot R model ODE predictions, Stan model predictions and data 
p3 <- odeall3 %>% ggplot2::ggplot() + 
  geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample)), shape = 5, size = 3, data = df4) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), shape = 20, data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p3)
#ggsave("featured.png",p3)


## ---- pred_test_stanv2 ----
#get predictions for first model  
x2 <- samp_m1 %>% select( contains("virus_pred2")) 
#observations of all samples for first individual only
x3 <- data.frame(x2[,1:Nobs[1]]) 
# turn into a long data frame that contains value for each individual and time
# some transformation shennanigans to get the values in the right order
df5 <- data.frame(time = rep(dat$time, nrow(samp_m1)), 
                  sample = rep(1:nrow(samp_m1), each = Nobs[1]),
                  value = as.vector(t(as.matrix(x3)))
)
#plot virus load curve for first individual from Stan model predictions
p4 <- ggplot2::ggplot() + 
  geom_point(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample)), shape = 5, size = 3, data = df4) +
  geom_point(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample)), shape = 15, size = 2, data = df5) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p4)



## ---- summarize_m1 ----
# get summaries across samples for all parameters
par_res <- res_m1$summary()


## ---- mean-par-test ----
# get means of 4 main parameters
# 24 values each (one for each individual) of each parameter
x <- par_res[1:(4*Nind),]$mean
odeall4 = NULL
for (n in 1:Nind){
  # virus load varies by dose, but barely so keeping it fixed here
  parms = c(alph=x[n],bet=x[n+Nind],gamm=x[n+2*Nind],et=x[n+2*Nind])
  # run model for each set of parameters for each individual
  oderes4 <- ode(y = y0, t = timevec, func = odemod1, parms = parms)
  # combine all results, also add individual ID as n and add dose for a given individual
  odeall4 = rbind(odeall4, cbind(oderes4, n, dat$dose_cat[match(n,dat$id)]))
}
odeall4 = data.frame(odeall4)

## ---- mean-pred-test ----
# linear model
vir1 = filter(par_res, grepl("virus_pred1",variable))
df6 = data.frame(y = vir1$mean, time = dat$time, n = dat$id, dose_cat = dat$dose_cat)
# log model
vir2 = filter(par_res, grepl("virus_pred2",variable))
df7 = data.frame(y = vir2$mean, time = dat$time, n = dat$id, dose_cat = dat$dose_cat)


#plot virus load curve for each individual from ODE model
p6 <- ggplot2::ggplot() + 
  geom_point(aes(x=time, y=y, group=as.factor(n),color=as.factor(n)), shape = 5, size = 3, data = df6) +
  geom_point(aes(x=time, y=y, group=as.factor(n),color=as.factor(n)), shape = 15, size = 2, data = df7) +
  geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n)), data = odeall4) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p6)
#ggsave("featured.png",p3)



## ---- make_stanmodel2 -----
# make Stan model. 
stanmod2 <- cmdstanr::cmdstan_model(stan_file = stanfile2, 
                                    compile = TRUE,
                                    pedantic=TRUE,
                                    force_recompile=TRUE)


## ---- fit-conditions ----
#settings for fitting
fs_m2 = list(warmup = 500,
             sampling = 500, 
             max_td = 18, #tree depth
             adapt_delta = 0.99,
             chains = 2,
             cores  = 2,
             seed = 1234,
             save_warmup = TRUE)

## ---- initialconditions2 ----
# separate definition of initial values
# use same as above, 
# but since/if we change number of chains we need to repeat this part
inits = NULL
for (n in 1:fs_m2$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m2$init = inits

## ---- priors-m2 --------
# values for prior distributions
pv2 = list(a0_mu = 1, a0_sd = 0.5,
          b0_mu = -1, b0_sd = 0.5,
          g0_mu = 0, g0_sd = 2,
          e0_mu = 1, e0_sd = 2
)


## ---- run_m2 ----
# adding priors to data
fitdatstan2 = c(fitdat,
               pv2,
               as = 12, #scaling of a0
               bs = 35, #scaling of b0
               modeltype = 1
)

res_m2 <- stanmod2$sample(data = fitdatstan2,
                          chains = fs_m2$chains,
                          init = fs_m2$init,
                          seed = fs_m2$seed,
                          parallel_chains = fs_m2$chains,
                          iter_warmup = fs_m2$warmup,
                          iter_sampling = fs_m2$sampling,
                          save_warmup = fs_m2$save_warmup,
                          max_treedepth = fs_m2$max_td,
                          adapt_delta = fs_m2$adapt_delta,
                          output_dir = filepath
)

## ---- savefits_m2 ----
res_m2$save_object(fs::path(filepath,filename))


## ---- run_m2_v2 ----
# adding priors to data
fitdatstan2 = c(fitdat,
                pv,
                as = 12, #scaling of a0
                bs = 35, #scaling of b0
                modeltype = 2
)

res_m2_v2 <- stanmod2$sample(data = fitdatstan2,
                          chains = fs_m2$chains,
                          init = fs_m2$init,
                          seed = fs_m2$seed,
                          parallel_chains = fs_m2$chains,
                          iter_warmup = fs_m2$warmup,
                          iter_sampling = fs_m2$sampling,
                          save_warmup = fs_m2$save_warmup,
                          max_treedepth = fs_m2$max_td,
                          adapt_delta = fs_m2$adapt_delta,
                          output_dir = filepath
)



## ---- savefits_m2_v2 ----
res_m2_v2$save_object(fs::path(filepath,filename))

## ---- loadfits --------
# loading previously saved fit.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local cloud-synced folder
# adjust accordingly for your setup
res_m2 <- readRDS(fs::path(filepath,filename))
res_m2_v2 <- readRDS(fs::path(filepath,filename))

## ---- diagnose_m2 ----
res_m2$cmdstan_diagnose()
res_m2_v2$cmdstan_diagnose()


## ---- get_samples_m2----
#this uses the posterior package to get draws
samp_m2 <- res_m2$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m2 <- res_m2$draws(inc_warmup = TRUE, format = "draws_df")



## ---- plot_par_m2 ----
# only a few parameters
plotpars = c("a0[1]","b0[1]","g0[1]","e0[1]","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(allsamp_m2, pars = plotpars)
bp3 <- bayesplot::mcmc_dens_overlay(samp_m2, pars = plotpars)
plot(bp1)
plot(bp3)


## ---- prep_data_m2 ----
# data manipulation to get in shape for plotting
# start with manipulation of posterior parameters
postdf1 <- samp_m2 %>% 
  select(!ends_with('prior')) %>% 
  select(contains("sigma")) 
# akward way of getting some further parameters
# namely values from first individual for a0,b0,g0,e0
postdf2 <- samp_m2 %>%
  select(contains("0[1]")) %>%
  rename_with(~ gsub("[1]", "", .x, fixed = TRUE) )
postdf <- cbind(postdf1, postdf2) 
postlong <- tidyr::pivot_longer(data = postdf, 
                                cols = everything() , 
                                names_to = "parname", 
                                values_to = "value") %>% 
  dplyr::mutate(type = "posterior")
# manipulation of prior parameters
priordf <-  samp_m2 %>% 
  select(ends_with('prior')) %>% 
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE) ) 
priorlong <- tidyr::pivot_longer(data = priordf, 
                                 cols = everything() , 
                                 names_to = "parname", 
                                 values_to = "value") %>% 
  dplyr::mutate(type = "prior")

ppdf <- dplyr::bind_rows(postlong,priorlong)

## ---- prior_post_m2 ----
m2_p1 <- ppdf %>%
  ggplot() +
  geom_density(aes(x = value, color = type), linewidth = 1) +
  facet_wrap("parname", scales = "free") +
  theme_minimal()
plot(m2_p1)


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
                      Estimate = mu2[,2],
                      Qmulo = mu2[,1], Qmuhi = mu2[,3],
                      Qsimlo = preds2[,1], Qsimhi = preds2[,3]
)

#make the plot
#not including the CI or prediction intervals since it looks too messy
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


