# code for fit exploration of tutorial part 4
# https://www.andreashandel.com/posts/longitudinal-multilevel-bayesian-analysis-4/


## ---- packages2 --------
library('dplyr') # for data manipulation
library('ggplot2') # for plotting
library('cmdstanr') #for model fitting
library('rethinking') #for model fitting
library('fs') #for file path


## ---- predict_function --------
# defining a function so I don't need to re-type the same code
predfunction <- function(fl,fitdat)
{

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
    #make variables for all versions of dose we use
    preddat = data.frame( id = sort(rep(seq(1,Ntot),length(timevec))),
                          time = rep(timevec,Ntot),
                          dose = 0,
                          dose_adj = 0,
                          dose_adj2 = 0,
                          dose_cat = "0",
                          dose_cat2 = 0
    )
    #add right dose information for each individual
    #this could likely be coded better
    for (k in 1:Ntot)
    {
      #get actual dose for a given individual
      #assign that dose
      #if statements because not every data set/model has each dose type
      if ("dose" %in% names(fitdat)) {
        nowdose = unique(fitdat$dose[fitdat$id == k])
        preddat[(preddat$id == k),"dose"] = nowdose
      }
      if ("dose_adj" %in% names(fitdat)) {
        nowdose_adj = unique(fitdat$dose_adj[fitdat$id == k])
        preddat[(preddat$id == k),"dose_adj"] = nowdose_adj
      }
      if ("dose_adj2" %in% names(fitdat)) {
        nowdose_adj2 = unique(fitdat$dose_adj2[fitdat$id == k])
        preddat[(preddat$id == k),"dose_adj2"] = nowdose_adj2
      }
      if ("dose_cat" %in% names(fitdat)) {
        nowdose_cat = unique(fitdat$dose_cat[fitdat$id == k])
        preddat[(preddat$id == k),"dose_cat"] = as.character(nowdose_cat)
      }
      if ("dose_cat2" %in% names(fitdat)) {
        nowdose_cat2 = unique(fitdat$dose_cat2[fitdat$id == k])
        preddat[(preddat$id == k),"dose_cat2"] = nowdose_cat2
      }
    }

    # pull out posterior samples for the parameters
    post <- extract.samples(nowmodel)

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
                              dose = preddat$dose_cat,
                              predtime = preddat$time,
                              Estimate = modmean,
                              Q79lo = modPI79[1,], Q79hi = modPI79[2,],
                              Q89lo = modPI89[1,], Q89hi = modPI89[2,],
                              Q97lo = modPI97[1,], Q97hi = modPI97[2,],
                              Qsimlo=modPIsim[1,], Qsimhi=modPIsim[2,]
    )
  } #end loop over all models
  return(fitpred)
} #end function computing predictions


## ---- plot_function --------
# defining a function so I don't need to re-type the same code
plotfunction <- function(fl,fitpred,fitdat)
{
  #list for storing all plots
  plotlist = vector(mode = "list", length = length(fl))

  #small data adjustment for plotting
  plotdat <- fitdat %>% data.frame()  %>%
    mutate(id = as.factor(id))  %>%
    mutate(dose = dose_cat)


  #looping over all models, creating and storing a plot for each
  for (n in 1:length(fl))
  {
    #adding titles to plots
    title = fl[[n]]$model

    plotlist[[n]] <- ggplot(data = fitpred[[n]], aes(x = predtime, y = Estimate, group = id, color = dose )) +
      geom_line(show.legend = F ) +
      geom_ribbon(aes(x=predtime, ymin=Q89lo, ymax=Q89hi, fill = dose, color = NULL), alpha=0.3, show.legend = F) +
      geom_ribbon(aes(x=predtime, ymin=Qsimlo, ymax=Qsimhi, fill = dose, color = NULL), alpha=0.1, show.legend = F) +
      geom_point(data = plotdat, aes(x = time, y = outcome, group = id, color = dose), shape = 1, size = 2, stroke = 2) +
      scale_y_continuous(limits = c(-30,50)) +
      labs(y = "Virus load",
           x = "days post infection") +
      theme_minimal() +
      ggtitle(title)
  }
  return(plotlist)
} #end function making plots


## ---- load_dat2 --------
#loading previously saved fits.
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_dat2", ext="Rds")
if (!file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits_dat2", ext="Rds")
}

fl <- readRDS(filepath)
fitdat <- fl[[1]]$fit@data

## ---- explore_dat2 --------
#Model 2a
print(precis(fl[[1]]$fit,depth=1),digits = 2)
#Model 4
a0mean = mean(precis(fl[[2]]$fit,depth=2,"a0")$mean)
b0mean = mean(precis(fl[[2]]$fit,depth=2,"b0")$mean)
print(precis(fl[[2]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))


## ---- compare_dat2 --------
comp <- compare(fl[[1]]$fit,fl[[2]]$fit)
print(comp)


## ---- predict_dat2 --------
fitpred <- predfunction(fl,fitdat)

## ---- plot_dat2 --------
plotlist <- plotfunction(fl,fitpred,fitdat)
plot(plotlist[[1]])
plot(plotlist[[2]])


## ---- load_big --------
#loading previously saved fits.
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_big", ext="Rds")
if (!file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits_big", ext="Rds")
}

fl <- readRDS(filepath)
fitdat <- fl[[1]]$fit@data

## ---- explore_big --------
a0mean = mean(precis(fl[[1]]$fit,depth=2,"a0")$mean)
b0mean = mean(precis(fl[[1]]$fit,depth=2,"b0")$mean)
print(precis(fl[[1]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))


## ---- predict_big --------
fitpred <- predfunction(fl,fitdat)

## ---- plot_big --------
plotlist <- plotfunction(fl,fitpred,fitdat)
#update plot
p1 <- plotlist[[1]] + scale_y_continuous(limits = c(-30,70))
plot(p1)
#save plot so we can use it in the blog post
ggsave(file = paste0("featured.png"), p1, dpi = 300, units = "in", width = 6, height = 6)



## ---- load_altpos --------
#loading previously saved fits.
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_altpos", ext="Rds")
if (!file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits_altpos", ext="Rds")
}
fl <- readRDS(filepath)
fitdat <- fl[[1]]$fit@data

## ---- explore_altpos_m4a --------
#Model 4a
a0mean = mean(precis(fl[[1]]$fit,depth=2,"a0")$mean)
b0mean = mean(precis(fl[[1]]$fit,depth=2,"b0")$mean)
print(precis(fl[[1]]$fit,depth=1),digits = 2)
print(c(a0mean,b0mean))

## ---- explore_altpos_m5 --------
#Model 5
print(precis(fl[[2]]$fit,depth=1),digits = 2)
a0mean = mean(precis(fl[[2]]$fit,depth=2,"a0")$mean)
b0mean = mean(precis(fl[[2]]$fit,depth=2,"b0")$mean)
print(c(a0mean,b0mean))

## ---- convert_altpos_m5 --------
# computing values that correspond to a1 and b1
a1est = (precis(fl[[2]]$fit,pars="a2")[1,]-1)*a0mean/max(fitdat$dose)
b1est = (precis(fl[[2]]$fit,pars="b2")[1,]-1)*b0mean/max(fitdat$dose)
print(c(a1est,b1est))

## ---- compare_altpos --------
compare(fl[[1]]$fit,fl[[2]]$fit)


## ---- predict_altpos --------
fitpred <- predfunction(fl,fitdat)

## ---- plot_altpos --------
plotlist <- plotfunction(fl,fitpred,fitdat)
plot(plotlist[[1]])
plot(plotlist[[2]])


## ---- load_cat --------
#loading previously saved fits.
filepath = fs::path("D:","Dropbox","datafiles","longitudinalbayes","ulamfits_cat", ext="Rds")
if (!file_exists(filepath))
{
  filepath = fs::path("C:","Data","Dropbox","datafiles","longitudinalbayes","ulamfits_cat", ext="Rds")
}
fl <- readRDS(filepath)
fitdat <- fl[[1]]$fit@data

## ---- explore_cat --------
#Model 4
print(precis(fl[[1]]$fit,depth=2),digits = 2)
#Model 6
print(precis(fl[[2]]$fit,depth=2),digits = 2)

## ---- compare_cat --------
compare(fl[[1]]$fit,fl[[2]]$fit)


## ---- predict_cat --------
fitpred <- predfunction(fl,fitdat)


## ---- plot_cat --------
plotlist <- plotfunction(fl,fitpred,fitdat)
plot(plotlist[[1]])
plot(plotlist[[2]])
