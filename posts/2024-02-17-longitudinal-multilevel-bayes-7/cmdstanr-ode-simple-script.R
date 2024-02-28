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

## ---- setup --------
############################################
# some general definitons and setup stuff
############################################
#setting random number seed for reproducibility
rngseed = 1234
# File locations and names
# adjust as needed
#filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes")
filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes")
filename = "cmdstanr-ode-simple.Rds"
stanfile <- here('posts','2024-02-17-longitudinal-multilevel-bayes-7',"stancode-ode-simple.stan")
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')




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
# loading data, path set in setup
simdat <- readRDS(simdatloc)
# using dataset 3 for fitting
# need to reformat all data as one list, this is how Stan needs it
# we used to have things arranged by observations at a given time for every individual
# for this model to make it work in Stan, we need to reorganize to
# have all data for one individual, then the next, etc.
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
fitdat=list(id=newdat$id, #an ID for each individual, for indexing
            outcome = newdat$outcome, #all outcomes
            time = newdat$time, #all times
            dose_adj = newdat$dose_adj,
            dose_level = as.numeric(factor(newdat$dose_adj)), 
            Ntot =  Ntot,
            Nobs =  Nobs,
            Nind = Nind,
            Ndose = Ndose,
            start = start,
            stop = stop,
            tstart = 0
)



# values for prior distributions
# allows for exploring different values without having to edit Stan model code
priorvals = list(a0_mu = 3, a0_sd = 1,
                 b0_mu = -1, b0_sd = 1,
                 g0_mu = 1, g0_sd = 1,
                 e0_mu = 1, e0_sd = 1,
                 V0_mu = exp(5), V0_sd = 10
)
fitdat = c(fitdat,priorvals)


## ---- make_stanmodel -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(stan_file = stanfile, 
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
set.seed(rngseed) #make inits reproducible
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

## ---- get_samples_m1 ----
#this uses the posterior package to get draws
samp_m1 <- res_m1$draws(inc_warmup = FALSE, format = "draws_df")
allsamp_m1 <- res_m1$draws(inc_warmup = TRUE, format = "draws_df")

## ---- diagnose_m1 ----
res_m1$cmdstan_diagnose()

## ---- summarize_m1 ----
# get summaries across samples for all parameters
par_res <- res_m1$summary()
print(head(par_res,15))


## ---- test_m1 ----
#get means of 4 main parameters
x <- par_res[3:(3+4*Nind),]$mean
t = seq(0,10,length=100)
odeall = NULL
for (n in 1:Nind){
  # virus load varies by dose, but barely so keeping it fixed here
  y0 = log(c(1e8,1,149))
  parms = exp(c(alph=x2[n],bet=x2[n+Nind],gamm=x2[n+2*Nind],et=x2[n+2*Nind]))
  # run model for each set of parameters for each individual
  oderes <- ode(y = y0, t = t, func = odemod, parms = parms)
  # combine all results, also add individual ID as n
  odeall=rbind(odeall, cbind(oderes,n))
}
#plot virus load curve for each individual
odeall = data.frame(odeall)
p1 <- odeall %>% ggplot2::ggplot() + geom_line(aes(x=time, y=X3, group=as.factor(n),color=as.factor(n)))
plot(p1)

x2 <- par_res %>% dplyr::filter(grepl("virus",variable))

## ---- plot_par_m1 ----
# only a few parameters
plotpars = c("a0[1]","b0[1]","g0[1]","e0[1]","sigma")
bayesplot::color_scheme_set("viridis")
bp1 <- bayesplot::mcmc_trace(samp_m1, pars = plotpars)
bp3 <- bayesplot::mcmc_dens_overlay(samp_m1, pars = plotpars)
plot(bp1)
plot(bp3)


## ---- prep_data_m1 ----
# data manipulation to get in shape for plotting
# data manipulation to get in shape for plotting
# start with manipulation of posterior parameters
postdf1 <- samp_m1 %>% 
  select(!ends_with('prior')) %>% 
  select(contains("sigma")) 
# akward way of getting some further parameters
# namely values from first individual for a0,b0,g0,e0
postdf2 <- samp_m1 %>%
  select(contains("0[1]")) %>%
  rename_with(~ gsub("[1]", "", .x, fixed = TRUE) )
postdf <- cbind(postdf1, postdf2) 
postlong <- tidyr::pivot_longer(data = postdf, 
                                cols = everything() , 
                                names_to = "parname", 
                                values_to = "value") %>% 
  dplyr::mutate(type = "posterior")
# manipulation of prior parameters
priordf <-  samp_m1 %>% 
  select(ends_with('prior')) %>% 
  rename_with(~ gsub("_prior", "", .x, fixed = TRUE) ) 
priorlong <- tidyr::pivot_longer(data = priordf, 
                                 cols = everything() , 
                                 names_to = "parname", 
                                 values_to = "value") %>% 
  dplyr::mutate(type = "prior")

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




