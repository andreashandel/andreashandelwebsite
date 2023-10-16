# Code only portion of tutorial part 1
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-1/
# This code simulates the data we'll fit in later posts

## ---- settings --------
## General settings
set.seed(123) #for reproducibility
# days at which we assume outcome is measured
timevec <- c(0.1,1,3,5,7,10,14,21,28,35,42)

#different number of individuals per dose to make it clearer which is which
#also, that's the structure of the data which motivated the tutorial
Nlow = 7; Nmed = 8; Nhigh = 9; filename = "simdat.Rds"
#if you want to explore how model fitting changes if you increase sample size
#turn on this line of code
#this is used in part 4 of the tutorial
#Nlow = 70; Nmed = 80; Nhigh = 90; filename = "simdat_big.Rds"

Ntot = Nlow + Nmed + Nhigh; #total number of individuals

# Set values for dose
# since we only consider dose on a log scale
# we'll log transform right here and then always use it in those log units
high_dose = log(1000)
med_dose = log(100)
low_dose = log(10)
dosevec = c(rep(low_dose,Nlow),rep(med_dose,Nmed),rep(high_dose,Nhigh))
# we are also creating a version of the dose variable
# that consists of ordered categories instead of numeric values
# we'll use that mostly for plotting
dosevec_cat = ordered(c(rep("low", Nlow),
                        rep("medium",Nmed),
                        rep("high",Nhigh)),
                      levels=c("low","medium","high"))


## Setting parameter values

## ---- commonpars --------
sigma = 1
a1 = 0.1
b1 = -0.1

## ---- m1pars --------
m1_mua = 3
m1_mub = 1
m1_sigmaa = 1
m1_sigmab = 1
m1_a0 = rnorm(n=Ntot, m1_mua, m1_sigmaa)
m1_b0 = rnorm(n=Ntot, m1_mub, m1_sigmaa)
# saving main parameters
m1pars = c(sigma = sigma, a1 = a1, b1 = b1,
           a0_mu = m1_mua, b0_mu = m1_mub)


## ---- m2pars --------
m2_mua = 3
m2_mub = 1
m2_sigmaa = 0.0001
m2_sigmab = 0.0001
m2_a0 = rnorm(n=Ntot, m2_mua, m2_sigmaa)
m2_b0 = rnorm(n=Ntot, m2_mub, m2_sigmab)
m2pars = c(sigma = sigma, a1 = a1, b1 = b1,
           a0_mu = m2_mua, b0_mu = m2_mub)


## ---- m3pars --------
m3_mua = 3
m3_mub = 1
m3_sigmaa = 0.1
m3_sigmab = 0.1
m3_a0 = rnorm(n=Ntot, m3_mua, m3_sigmaa)
m3_b0 = rnorm(n=Ntot, m3_mub, m3_sigmaa)
m3pars = c(sigma = sigma, a1 = a1, b1 = b1,
           a0_mu = m3_mua, b0_mu = m3_mub)


## ---- m1sims --------
m1_alpha = m1_a0 + a1*(dosevec - med_dose)
m1_beta = m1_b0 + b1*(dosevec - med_dose)
#doing matrix multiplication to get time-series for each individual
#for that to work, the timevec vector needs to be transposed
m1_mu =  exp(m1_alpha) %*% t(log(timevec)) - exp(m1_beta) %*% t(timevec)
# apply variation following a normal distribution to each data point
m1_y = rnorm(length(m1_mu),m1_mu, sigma)
# in a final step, we reorganize the data into a long data frame with
# columns id, time, dose, model,
# the deterministic mean mu, and the normally distributed outcome.
# We store dose in 3 versions, the original (log transformed one),
# the one that has the middle value subtracted, and a categorical one.
# Note that trick using sort to get time in the right order.
# Not a robust way of doing things, but works here
m1_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)),
                     mu = as.vector(m1_mu),
                     outcome = as.vector(m1_y),
                     model = "m1")


## ---- m2m3sims --------
#model 2
m2_alpha = m2_a0 + a1*(dosevec - med_dose)
m2_beta = m2_b0 + b1*(dosevec - med_dose)
m2_mu =  exp(m2_alpha) %*% t(log(timevec)) - exp(m2_beta) %*% t(timevec)
m2_y = rnorm(length(m2_mu),m2_mu, sigma)
m2_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)),
                     mu = as.vector(m2_mu),
                     outcome = as.vector(m2_y),
                     model = "m2")

#model 3
m3_alpha = m3_a0 + a1*(dosevec - med_dose)
m3_beta = m3_b0 + b1*(dosevec - med_dose)
m3_mu =  exp(m3_alpha) %*% t(log(timevec)) - exp(m3_beta) %*% t(timevec)
m3_y = rnorm(length(m3_mu),m3_mu, sigma)
m3_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)),
                     dose_adj = rep(dosevec,length(timevec))-med_dose,
                     dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)),
                     mu = as.vector(m3_mu),
                     outcome = as.vector(m3_y),
                     model = "m3")



## ---- packages --------
library('ggplot2')

## ---- makeplots --------
p1 <- ggplot(m1_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,200)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)") +
  theme_minimal()


p2 <- ggplot(m2_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)") +
  theme_minimal()

p3 <- ggplot(m3_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)") +
  theme_minimal()

## ---- showplots --------
plot(p1)
plot(p2)
plot(p3)

## ---- savesims --------
#save a plot so we can use it in the blog post
simdat <- list(m1 = m1_dat, m2 = m2_dat, m3 = m3_dat, m1pars = m1pars, m2pars = m2pars, m3pars = m3pars, Nlow = Nlow, Nmed = Nmed, Nhigh = Nhigh)
saveRDS(simdat, file = filename)
ggsave(file = paste0("featured.png"), p3, dpi = 300, units = "in", width = 6, height = 6)

