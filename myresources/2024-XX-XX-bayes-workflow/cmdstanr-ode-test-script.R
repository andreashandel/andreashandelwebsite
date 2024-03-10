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
#filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes")
filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes")
filename = "cmdstanr-ode-testing.Rds"
stanfile <- here('myresources','2024-XX-XX-bayes-workflow',"stancode-ode-testing.stan")
stanfile2 <- here('myresources','2024-XX-XX-bayes-workflow',"stancode-ode-testing2.stan")
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')

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
                 g0_mu = 0, g0_sd = 0.5,
                 e0_mu = 1, e0_sd = 0.5
)


## ---- r-simulation --------
Nsamp = 20 #number of samples
timevec = seq(0,max(dat$time),length=1000) #time
odeall1 = NULL #set up empty array to hold all simulations
odeall2 = NULL #same empty array for model version 2
U0 = 1e12; #number of uninfected cells
I0 = 1e-15; #set to small non-zero so we can take log
V0 = 1e-3; #same value for everyone - doesn't agree with data but easy to try
# scaling such that parameter distributions are close to 1
as = 12;
bs = 35;
R0 = rep(0,Nsamp) #record R0 for all samples
set.seed(rngseed) #set seed for reproducibility
# could be done more efficiently without a loop
for (n in 1:Nsamp)
{
  alph = exp(as*rnorm(1,pv$a0_mu, pv$a0_sd))
  bet = exp(bs*rnorm(1,pv$b0_mu, pv$b0_sd))
  gamm = exp(rnorm(1,pv$g0_mu, pv$g0_sd)) 
  et = exp(rnorm(1,pv$e0_mu, pv$e0_sd))
  R0[n] = bet*alph*U0/(gamm*et)
  odeparms = c(alph=alph,bet=bet,gamm=gamm,et=et)  
  # Uninfected cells, infected cells, virus
  # note that I'm starting with 1 infected cells so I can take the log
  y0 = c(U0,I0,V0)  
  # run the ODE model twice, in linear and log space
  # should give the same results 
  oderes1 <- ode(y = y0, t = timevec, func = odemod1, parms = odeparms, method = "bdf", atol = 1e-10, rtol = 1e-10)
  oderes2 <- ode(y = log(y0), t = timevec, func = odemod2, parms = odeparms, method = "bdf")
  # combine all results, also add individual ID as n and add dose for a given individual
  odeall1=rbind(odeall1, cbind(oderes1,n))
  odeall2=rbind(odeall2, cbind(oderes2,n))
}
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





## ---- make_stanmodel1 -----
# make Stan model. 
stanmod1 <- cmdstanr::cmdstan_model(stan_file = stanfile, 
                                    compile = TRUE,
                                    pedantic=TRUE,
                                    force_recompile=TRUE)

## ---- show_stancode -----
print(stanmod1)

## ---- fitconditions ----
#settings for fitting
fs_m1 = list(warmup = 1,
             sampling = 40, 
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
init_vals_1chain <- function() (list(a0 = rep(pv$a0_mu,Nind),
                                     b0 = rep(pv$b0_mu,Nind),
                                     g0 = rep(pv$g0_mu,Nind),
                                     e0 = rep(pv$e0_mu,Nind),
                                     sigma = 1
))
#set.seed(rngseed) #make inits reproducible
# init_vals_1chain <- function() (list(a0 = runif(Nind,5,5),
#                                      b0 = runif(Nind,-22,-22),
#                                      g0 = runif(Nind,-1,-1),
#                                      e0 = runif(Nind,0.4,0.5),
#                                      sigma = runif(1,0,1)
#                                      ))

inits = NULL
for (n in 1:fs_m1$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m1$init = inits


## ---- run_m1 ----
fitdatstan = c(fitdat,pv)
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


## ---- par_test_ind1 ----
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
p3 <- odeall3 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p3)
#ggsave("featured.png",p3)


## ---- pred_test_ind1 ----
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
#plot virus load curve for first individual from Stan model predictions
p4 <- df4 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p4)

## ---- pred_test_ind1v2 ----
#get predictions for first model  
x2 <- samp_m1 %>% select( contains("virus_pred2")) 
#observations of all samples for first individual only
x3 <- data.frame(x2[,1:Nobs[1]]) 
# turn into a long data frame that contains value for each individual and time
# some transformation shennanigans to get the values in the right order
df4 <- data.frame(time = rep(dat$time, nrow(samp_m1)), 
                  sample = rep(1:nrow(samp_m1), each = Nobs[1]),
                  value = as.vector(t(as.matrix(x3)))
)
#plot virus load curve for first individual from Stan model predictions
p5 <- df4 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p5)


## ---- summarize_m1 ----
# get summaries across samples for all parameters
par_res <- res_m1$summary()


## ---- par_test_m1 ----
#get means of 4 main parameters
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
#plot virus load curve for each individual from ODE model
odeall4 = data.frame(odeall4)
p6 <- odeall4 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme(legend.position="none")
plot(p6)
#ggsave("featured.png",p3)


## ---- pred_test_m1 ----
vir1 = filter(par_res, grepl("virus_pred1",variable))
df1 = data.frame(y = vir1$mean, x = dat$time, n = dat$id, dose_cat = dat$dose_cat)
p7 <- df1 %>% 
      ggplot2::ggplot() +
      geom_line(aes(x=x,y=y, group=as.factor(n),color=dose_cat)) +
      geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
      scale_y_continuous(limits = c(-30,50)) + 
      scale_x_continuous(limits = c(0,45)) +
      theme(legend.position="none")
plot(p7)

vir2 = filter(par_res, grepl("virus_pred2",variable))
df2 = data.frame(y = vir2$mean, x = dat$time, n = dat$id, dose_cat = dat$dose_cat)
p5 <- df2 %>% 
  ggplot2::ggplot() + 
  geom_line(aes(x=x,y=y, group=as.factor(n),color=dose_cat)) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme(legend.position="none")
plot(p5)

## ---- make_stanmodel2 -----
# make Stan model. 
stanmod2 <- cmdstanr::cmdstan_model(stan_file = stanfile2, 
                                    compile = TRUE,
                                    pedantic=TRUE,
                                    force_recompile=TRUE)


## ---- fitconditions ----
#settings for fitting
fs_m2 = list(warmup = 500,
             sampling = 500, 
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
inits = NULL
for (n in 1:fs_m2$chains)
{
  inits[[n]] = init_vals_1chain()
}
fs_m2$init = inits


## ---- run_m2 ----
res_m2 <- stanmod2$sample(data = fitdatstan,
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


## ---- savefits ----
res_m2$save_object(fs::path(filepath,filename))

## ---- loadfits --------
# loading previously saved fit.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local cloud-synced folder
# adjust accordingly for your setup
res_m2 <- readRDS(fs::path(filepath,filename))


## ---- get_samples_m2----
#this uses the posterior package to get draws
samp_m2 <- res_m2$draws(inc_warmup = FALSE, format = "draws_df")


## ---- par_test_ind2 ----
#get samples of parameters for first individual 
x1 <- samp_m2 %>% select(c("alph[1]","bet[1]","gamm[1]","et[1]")) 
odeall3 = NULL
for (n in 1:nrow(x1)){
  # virus load varies by dose, but barely so keeping it fixed here
  y0 = c(U0,1,V0)
  parms = c(alph=x1$`alph[1]`[n],bet=x1$`bet[1]`[n],gamm=x1$`gamm[1]`[n],et=x1$`et[1]`[n])
  # run model for each set of parameters for each individual
  oderes3 <- ode(y = y0, t = timevec, func = odemod1, parms = parms)
  # combine all results, also add individual ID as n and add dose for a given individual
  odeall3 = rbind(odeall3, cbind(oderes3, n, dat$dose_cat[match(n,dat$id)]))
}
#plot virus load curve for each individual from ODE model
odeall3 = data.frame(odeall3)
p3 <- odeall3 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p3)
#ggsave("featured.png",p3)


## ---- pred_test_ind1 ----
#get predictions for first model  
x2 <- samp_m2 %>% select( contains("virus_pred")) 
#observations of all samples for first individual only
x3 <- data.frame(x2[,1:Nobs[1]]) 
# turn into a long data frame that contains value for each individual and time
# some transformation shennanigans to get the values in the right order
df4 <- data.frame(time = rep(dat$time, nrow(samp_m2)), 
                  sample = rep(1:nrow(samp_m2), each = Nobs[1]),
                  value = as.vector(t(as.matrix(x3)))
)
#plot virus load curve for first individual from Stan model predictions
p4 <- df4 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=value, group=as.factor(sample),color=as.factor(sample))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme_minimal() +
  theme(legend.position="none")
plot(p4)


## ---- summarize_m2 ----
# get summaries across samples for all parameters
par_res <- res_m2$summary()


## ---- par_test_m2 ----
#get means of 4 main parameters
# 24 values each (one for each individual) of each parameter
x <- par_res[1:(4*Nind),]$mean
odeall4 = NULL
for (n in 1:Nind){
  # virus load varies by dose, but barely so keeping it fixed here
  y0 = c(U0,1,V0)
  parms = c(alph=x[n],bet=x[n+Nind],gamm=x[n+2*Nind],et=x[n+2*Nind])
  # run model for each set of parameters for each individual
  oderes4 <- ode(y = y0, t = timevec, func = odemod1, parms = parms)
  # combine all results, also add individual ID as n and add dose for a given individual
  odeall4 = rbind(odeall4, cbind(oderes4, n, dat$dose_cat[match(n,dat$id)]))
}
#plot virus load curve for each individual from ODE model
odeall4 = data.frame(odeall4)
p6 <- odeall4 %>% ggplot2::ggplot() + geom_line(aes(x=time, y=log(X3), group=as.factor(n),color=as.factor(n))) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme(legend.position="none")
plot(p6)
#ggsave("featured.png",p3)


## ---- pred_test_m1 ----
vir1 = filter(par_res, grepl("virus_pred1",variable))
df1 = data.frame(y = vir1$mean, x = dat$time, n = dat$id, dose_cat = dat$dose_cat)
p7 <- df1 %>% 
  ggplot2::ggplot() +
  geom_line(aes(x=x,y=y, group=as.factor(n),color=dose_cat)) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme(legend.position="none")
plot(p7)

vir2 = filter(par_res, grepl("virus_pred2",variable))
df2 = data.frame(y = vir2$mean, x = dat$time, n = dat$id, dose_cat = dat$dose_cat)
p5 <- df2 %>% 
  ggplot2::ggplot() + 
  geom_line(aes(x=x,y=y, group=as.factor(n),color=dose_cat)) +
  geom_point(aes(x=time,y=outcome,group=as.factor(id),color=as.factor(dose_cat)), data = dat) +
  scale_y_continuous(limits = c(-30,50)) + 
  scale_x_continuous(limits = c(0,45)) +
  theme(legend.position="none")
plot(p5)

