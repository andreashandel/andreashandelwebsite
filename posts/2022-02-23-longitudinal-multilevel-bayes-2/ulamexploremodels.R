# part 2 of code of ulam/rethinking tutorial (part 2)
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-2/
# this has the code part that explores the fits produced by the previous R script


## ---- packages2 --------
library('cmdstanr') #for model fitting
library('fs') #for file path
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('rstan') #is apparently called by some other functions in cmdstanr
library('rethinking') #for model fitting

## ---- loadfits --------
# loading list of previously saved fits.
# useful if we don't want to re-fit
# every time we want to explore the results.
# since the file is too large for GitHub
# it is stored in a local folder
# adjust accordingly for your setup
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
if (!fs::file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits", ext="Rds")
}

fl <- readRDS(filepath)
# also load data file used for fitting
simdatloc <- here::here('posts','2022-02-22-longitudinal-multilevel-bayes-1','simdat.Rds')
simdat <- readRDS(simdatloc)
#pull our the data set we used for fitting
#if you fit a different one of the simulated datasets, change accordingly
fitdat <- simdat$m3
#contains parameters used for fitting
pars <- simdat$m3pars

## ---- diagnostics ------
# Model 2a summary
show(fl[[5]]$fit)

## ---- traceplot ------
# Model 2a trace plots
# for some reason didn't work on last compile
rethinking::traceplot(fl[[5]]$fit, pars = c("a0","b0","a1","b1","sigma"))

## ---- trankplot ------
# Model 2a trank plots
rethinking::trankplot(fl[[5]]$fit, pars = c("a0","b0","a1","b1","sigma"))

## ---- pairplot ------
# Model 2a pair plot
# Correlation between posterior samples of parameters
rethinking::pairs(fl[[5]]$fit, pars = c("a0","b0","a1","b1","sigma"))


## ---- mod_1_3_prior --------
#get priors and posteriors for models 1 and 3
m1prior <- rethinking::extract.prior(fl[[1]]$fit, n = 1e4)
m1post <- rethinking::extract.samples(fl[[1]]$fit, n = 1e4)

m3prior <- rethinking::extract.prior(fl[[3]]$fit, n = 1e4)
m3post <- rethinking::extract.samples(fl[[3]]$fit, n = 1e4)

## ---- mod_1_3_prior_plots --------
#showing density plots for a0
plot(density(m1prior$a0), xlim = c (-20,20), ylim = c(0,2), lty=2)
lines(density(m1post$a0), lty=1)
lines(density(m3prior$a0), col = "blue", lty=2)
lines(density(m3post$a0), col = "blue", lty=1)

#showing density plots for b0
plot(density(m1prior$b0), xlim = c (-20,20), ylim = c(0,2), lty=2)
lines(density(m1post$b0), lty=1)
lines(density(m3prior$b0), col = "blue", lty=2)
lines(density(m3post$b0), col = "blue", lty=1)

#showing density plots for a1
plot(density(m1prior$a1), xlim = c (-3,3), ylim = c(0,2), lty=2)
lines(density(m1post$a1), lty=1)
lines(density(m3prior$a1), col = "blue", lty=2)
lines(density(m3post$a1), col = "blue", lty=1)

#showing density plots for b1
plot(density(m1prior$b1), xlim = c (-3,3), ylim = c(0,2), lty=2)
lines(density(m1post$b1), lty=1)
lines(density(m3prior$b1), col = "blue", lty=2)
lines(density(m3post$b1), col = "blue", lty=1)


## ---- mod_1_3_pair_plots --------
# all "a" parameters - too big to show
#pairs(fl[[1]]$fit, pars = c("a0","a1"))
# a few parameters for each dose
#low dose
rethinking::pairs(fl[[1]]$fit, pars = c("a0[1]","a0[2]","a0[3]","a0[4]","a1"))
#medium dose
rethinking::pairs(fl[[1]]$fit, pars = c("a0[8]","a0[9]","a0[10]","a0[11]","a1"))
#high dose
rethinking::pairs(fl[[1]]$fit, pars = c("a0[16]","a0[17]","a0[18]","a0[19]","a1"))



## ---- mod_1_3_exploration --------
# Model 1
a0mean = mean(rethinking::precis(fl[[1]]$fit,depth=2,"a0")$mean)
b0mean = mean(rethinking::precis(fl[[1]]$fit,depth=2,"b0")$mean)
print(rethinking::precis(fl[[1]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))

# Model 3
a0mean = mean(rethinking::precis(fl[[3]]$fit,depth=2,"a0")$mean)
b0mean = mean(rethinking::precis(fl[[3]]$fit,depth=2,"b0")$mean)
print(rethinking::precis(fl[[3]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))




## ---- mod_2_2a_exploration -------
# Compare models 2 and 2a
# first we compute the mean across individuals for model 2
a0mean = mean(rethinking::precis(fl[[2]]$fit,depth=2,"a0")$mean)
b0mean = mean(rethinking::precis(fl[[2]]$fit,depth=2,"b0")$mean)

#rest of model 2
print(rethinking::precis(fl[[2]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))
#model 2a
print(rethinking::precis(fl[[5]]$fit,depth=1),digits = 2)


## ---- mod_comparison --------
rethinking::compare(fl[[1]]$fit,fl[[3]]$fit,fl[[2]]$fit,fl[[5]]$fit)


## ---- mod_4_4a_prior --------
#get priors and posteriors for models 4 and 4a
m4prior <- rethinking::extract.prior(fl[[4]]$fit, n = 1e4)
m4post <- rethinking::extract.samples(fl[[4]]$fit, n = 1e4)

m4aprior <- rethinking::extract.prior(fl[[6]]$fit, n = 1e4)
m4apost <- rethinking::extract.samples(fl[[6]]$fit, n = 1e4)

## ---- mod_4_4a_prior_plots --------
#showing density plots for a0
plot(density(m4prior$mu_a), xlim = c (-10,10), ylim = c(0,2), lty=2)
lines(density(m4post$mu_a), lty=1)
lines(density(m4aprior$mu_a), col = "blue", lty=2)
lines(density(m4apost$mu_a), col = "blue", lty=1)

#showing density plots for b0
plot(density(m4prior$mu_b), xlim = c (-10,10), ylim = c(0,2), lty=2)
lines(density(m4post$mu_b), lty=1)
lines(density(m4aprior$mu_b), col = "blue", lty=2)
lines(density(m4apost$mu_b), col = "blue", lty=1)

#showing density plots for a1
plot(density(m4prior$a1), xlim = c (-3,3), ylim = c(0,2), lty=2)
lines(density(m4post$a1), lty=1)
lines(density(m4aprior$a1), col = "blue", lty=2)
lines(density(m4apost$a1), col = "blue", lty=1)

#showing density plots for b1
plot(density(m4prior$b1), xlim = c (-3,3), ylim = c(0,2), lty=2)
lines(density(m4post$b1), lty=1)
lines(density(m4aprior$b1), col = "blue", lty=2)
lines(density(m4apost$b1), col = "blue", lty=1)


## ---- mod_4_4a_pair_plots --------
# a few parameters for each dose
#low dose
rethinking::pairs(fl[[4]]$fit, pars = c("a0[1]","a0[2]","a0[3]","a0[4]","a1"))
#medium dose
rethinking::pairs(fl[[4]]$fit, pars = c("a0[8]","a0[9]","a0[10]","a0[11]","a1"))
#high dose
rethinking::pairs(fl[[4]]$fit, pars = c("a0[16]","a0[17]","a0[18]","a0[19]","a1"))
# mean of a0 prior
rethinking::pairs(fl[[4]]$fit, pars = c("mu_a","mu_b","a1","b1"))

#saving one plot so I can use as featured image
png(filename = "featured.png", width = 6, height = 6, units = "in", res = 300)
rethinking::pairs(fl[[4]]$fit, pars = c("mu_a","mu_b","a1","b1"))
dev.off()


## ---- mod_4_4a_exploration --------
# model 4
print(rethinking::precis(fl[[4]]$fit,depth=1),digits = 2)
# model 4a
print(precis(fl[[6]]$fit,depth=1),digits = 2)


## ---- mod_4_4a_comparison --------
rethinking::compare(fl[[3]]$fit,fl[[4]]$fit,fl[[6]]$fit, func = WAIC)
rethinking::compare(fl[[3]]$fit,fl[[4]]$fit,fl[[6]]$fit, func = PSIS)





## ---- mod_5_pair_plots --------
# a few parameters for each dose
#low dose
rethinking::pairs(fl[[7]]$fit, pars = c("a0[1]","a0[2]","a0[3]","a0[4]","a0[5]"))
#medium dose
rethinking::pairs(fl[[7]]$fit, pars = c("a0[8]","a0[9]","a0[10]","a0[11]","a0[12]"))
#high dose
rethinking::pairs(fl[[7]]$fit, pars = c("a0[16]","a0[17]","a0[18]","a0[19]","a0[20]"))


## ---- mod_5_exploration --------
a0mean = mean(rethinking::precis(fl[[7]]$fit,depth=2,"a0")$mean)
b0mean = mean(rethinking::precis(fl[[7]]$fit,depth=2,"b0")$mean)
print(rethinking::precis(fl[[7]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))


## ---- mod_5_comparison --------
rethinking::compare(fl[[3]]$fit,fl[[4]]$fit,fl[[7]]$fit)




## ---- computepredictions --------
#small data adjustment for plotting
plotdat <- fitdat %>% data.frame()  %>%
                      mutate(id = as.factor(id))  %>%
                      mutate(dose = dose_cat)

#this will contain all the predictions from the different models
fitpred = vector(mode = "list", length = length(fl))

# we are looping over each fitted model
for (n in 1:length(fl))
{
  #get current model
  nowmodel = fl[[n]]$fit

  #make new data for which we want predictions
  #specifically, more time points so the curves are smoother
  timevec = seq(from = 0.1, to = max(fitdat$time), length=100)
  Ntot = max(fitdat$id)
  #new data used for predictions
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

  # pull out posterior samples for the parameters
  post <- rethinking::extract.samples(nowmodel)

  # estimate and CI for parameter variation
  # this uses the link function from rethinking
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
  #and store in a list for each model
  fitpred[[n]] = data.frame(id = as.factor(preddat$id),
                            dose = as.factor(preddat$dose_cat),
                            predtime = preddat$time,
                            Estimate = modmean,
                            Q79lo = modPI79[1,], Q79hi = modPI79[2,],
                            Q89lo = modPI89[1,], Q89hi = modPI89[2,],
                            Q97lo = modPI97[1,], Q97hi = modPI97[2,],
                            Qsimlo=modPIsim[1,], Qsimhi=modPIsim[2,]
                            )
} #end loop over all models



## ---- makeplots --------
#list for storing all plots
plotlist = vector(mode = "list", length = length(fl))

#looping over all models, creating and storing a plot for each
for (n in 1:length(fl))
{
  #adding titles to plots
  title = fl[[n]]$model

  plotlist[[n]] <- ggplot(data = fitpred[[n]], aes(x = predtime, y = Estimate, group = id, color = dose ) ) +
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
}



## ---- mod_1_3_plots --------
plot(plotlist[[1]])
plot(plotlist[[3]])


## ---- mod_2_2a_plots --------
plot(plotlist[[2]])
plot(plotlist[[5]])


## ---- mod_4_4a_plots --------
plot(plotlist[[4]])
plot(plotlist[[6]])


## ---- mod_5_plots --------
plot(plotlist[[7]])


## ---- additional-code -------
for (n in 1:length(fl)) {print(fl[[n]]$runtime)}

