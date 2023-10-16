# code only portion of tutorial part 1
# Simulating data

## ---- settings --------
## General settings
set.seed(123) #for reproducibility
# days at which we assume outcome is measured
timevec <- c(0.1,1,3,5,7,10,14,21,28,35,42)

#do different number of individuals per dose
#just to make it clearer which is which and also since that's the data structure
Nlow = 7; Nmed = 8; Nhigh = 9;
#if you want to explore how model fitting changes if you increase sample size
#turn on this line of code
#Nlow = 70; Nmed = 80; Nhigh = 90;

Ntot = Nlow + Nmed + Nhigh; #total number of individuals

# Set values for dose
# since we only consider dose on a log scale
# we'll log transform right here and then always use it in those log units
high_dose = log(1000)
med_dose = log(100)
low_dose = log(10)
dosevec = c(rep(low_dose,Nlow),rep(med_dose,Nmed),rep(high_dose,Nhigh))
# we are also creating a version that is just ordered categories instead of numeric values
# we'll use that for plotting and in the alternative analysis in a separate post
dosevec_cat = ordered(c(rep("low", Nlow),rep("medium",Nmed),rep("high",Nhigh)),levels=c("low","medium","high"))


## Setting parameter values

## ---- commonpars --------
sigma = 1
a1 = 0.2
b1 = -0.2

## ---- m1pars --------
m1_mua = 3
m1_mub = 1
m1_sigmaa = 1
m1_sigmab = 1
m1_a0 = rnorm(n=Ntot, m1_mua, m1_sigmaa)
m1_b0 = rnorm(n=Ntot, m1_mub, m1_sigmaa)
# saving main parameters
m1pars = c(sigma = sigma, a1 = a1, b1 = b1, a0_mu = m1_mua, b0_mu = m1_mub)

## ---- m2pars --------
m2_mua = 3
m2_mub = 1
m2_sigmaa = 0.0001
m2_sigmab = 0.0001
m2_a0 = rnorm(n=Ntot, m2_mua, m2_sigmaa)
m2_b0 = rnorm(n=Ntot, m2_mub, m2_sigmab)
m2pars = c(sigma = sigma, a1 = a1, b1 = b1, a0_mu = m2_mua, b0_mu = m2_mub)

## ---- m3pars --------
m3_mua = 3
m3_mub = 1
m3_sigmaa = 0.1
m3_sigmab = 0.1
m3_a0 = rnorm(n=Ntot, m3_mua, m3_sigmaa)
m3_b0 = rnorm(n=Ntot, m3_mub, m3_sigmaa)
m3pars = c(sigma = sigma, a1 = a1, b1 = b1, a0_mu = m3_mua, b0_mu = m3_mub)

## ---- m1sims --------
m1_alpha = m1_a0 + a1*(log(dosevec) - log(med_dose))
m1_beta = m1_b0 + b1*(log(dosevec) - log(med_dose))
#doing matrix multiplication to get time-series for each individual
#for that to work, the timevec vector needs to be transposed
m1_mu =  exp(m1_alpha) %*% t(log(timevec)) - exp(m1_beta) %*% t(timevec)
# apply variation following a normal distribution to each individual
m1_y = apply(m1_mu,2,rnorm,sigma)
# in a final step, we reorganize the data into a long data frame with
# columns id, time, dose, model,
# the deterministic mean mu, and the normally distributed outcome
# we store dose in 3 versions, the original (log transformed one), the one that has the middle value subtracted, and a categorical one
# note that little trick using sort to get time in the right order
# not a robust way of doing things, but works here
m1_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)), dose_adj = rep(dosevec,length(timevec))-med_dose, dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)), model = "m1",
                     mu = as.vector(m1_mu), outcome = rnorm(n=length(m1_mu),mean=as.vector(m1_mu),sd=sigma)  )


## ---- m2m3sims --------
#model 2
m2_alpha = m2_a0 + a1*(log(dosevec) - log(med_dose))
m2_beta = m2_b0 + b1*(log(dosevec) - log(med_dose))
m2_mu =  exp(m2_alpha) %*% t(log(timevec)) - exp(m2_beta) %*% t(timevec)
m2_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)), dose_adj = rep(dosevec,length(timevec))-med_dose, dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)), model = "m2",
                     mu = as.vector(m2_mu), outcome = rnorm(n=length(m2_mu),mean=as.vector(m2_mu),sd=sigma)  )

#model 3
m3_alpha = m3_a0 + a1*(log(dosevec) - log(med_dose))
m3_beta = m3_b0 + b1*(log(dosevec) - log(med_dose))
m3_mu =  exp(m3_alpha) %*% t(log(timevec)) - exp(m3_beta) %*% t(timevec)
m3_y = apply(m3_mu,2,rnorm,sigma)
m3_dat <- data.frame(id = rep(1:Ntot,length(timevec)),
                     dose = rep(dosevec,length(timevec)), dose_adj = rep(dosevec,length(timevec))-med_dose, dose_cat =  rep(dosevec_cat,length(timevec)),
                     time = sort(rep(timevec,Ntot)), model = "m3",
                     mu = as.vector(m3_mu), outcome = rnorm(n=length(m3_mu),mean=as.vector(m3_mu),sd=sigma)  )



## ---- packages --------
library('ggplot2')

## ---- makeplots --------
p1 <- ggplot(m1_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,150)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)")

p2 <- ggplot(m2_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)")

p3 <- ggplot(m3_dat) +
  geom_line(aes(x=time, y=mu, col = dose_cat, group = id)) +
  geom_point(aes(x=time, y=outcome, col = dose_cat)) +
  scale_y_continuous(limits = c(-30,50)) +
  labs(y = "Outcome (log virus load)",  x = "Time (days)")

## ---- showplots --------
plot(p1)
plot(p2)
plot(p3)

## ---- savesims --------
ggsave(file = paste0("simdata_m3.png"), p3, dpi = 300, units = "in", width = 7, height = 7)
simdat <- list(m1 = m1_dat, m2 = m2_dat, m3 = m3_dat, m1pars = m1pars, m2pars = m2pars, m3pars = m3pars)
saveRDS(simdat, "simdat.Rds")


