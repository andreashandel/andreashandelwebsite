# code that generates cholesterol data

## ---- packages --------
library('ggplot2')

## ---- settings --------
## General settings
set.seed(123) #for reproducibility
N = 500 #medium-sized dataset

## ---- datageneration --------
# patient ID, just a number
pid = 1:N
# dose comes in discrete units, we might treat it continuous
# include 0 to indicate some people are not treated
dose = sample(c(0,5,10,20,40),N,replace = TRUE)
# assume all individuals on statins are uniformly distributed between those ages
age = runif(N,min=30,max=100)
# BMI that seems somewhat close to population, log-normal because it needs to be positive
bmi = rlnorm(N,meanlog=log(25),sdlog=log(1.2))
# Sex, simply M/F assume not quite 50/50 distribution
sex = as.factor(ifelse(rbinom(N, size = 1, prob = 0.6),"M","F"))
# previous condition, code as Yes/No
previous = as.factor(ifelse(rbinom(N, size = 1, prob = 0.8),"No","Yes"))
# study site - an ID for 3 sites. Assume half at site 1, fewer at others
site = as.factor(sample(c("S1","S2","S3"),N, replace = TRUE, prob = c(0.5, 0.3, 0.2)))
# LDL cholesterol values
# assuming drug dose lowers ldl in a exponential/saturing manner, and ldl increases with squareroot of age
ldl = 100 - 10*(1-exp(-0.1*dose)) + 5*sqrt(age)
#stick it all into a data frame
cholesteroldata = data.frame(pid=pid,ldl=ldl,dose=dose,age=age,bmi=bmi,sex=sex,previous=previous,site=site)





## ---- showplots --------
plot(cholesteroldata$dose,cholesteroldata$ldl)
plot(cholesteroldata$age,cholesteroldata$ldl)
plot(cholesteroldata$bmi,cholesteroldata$ldl)
boxplot(ldl ~ sex, data = cholesteroldata)
boxplot(ldl ~ site, data = cholesteroldata)
boxplot(ldl ~ previous, data = cholesteroldata)



## ---- savesims --------
write.csv(cholesteroldata,file = here::here("data/cholesteroldata.csv"),row.names = FALSE)


