########################################################
# Model fitting with Turing/Julia 
# Fitting longitudinal data with a Bayesian hierarchical model
# Last updated 2024-02-16 (or later) by Andreas Handel
########################################################

######################################
# set directory
######################################
cd(@__DIR__)
pwd() #get/check directory

## ---- loadpackages --------
######################################
# install/load all needed packages
# https://docs.julialang.org/en/v1/stdlib/Pkg/
######################################
using Pkg; # load the package manager package
Pkg.activate(".") #activate project and package manager
using Random #for random numbers
using Distributions #for probabilistic distributions
using Turing #Bayesian fitting routines
using Chain #for chaining commands like the R pipe
using MCMCChains #for visualization of turing output
using Plots # for basic plots
using StatsPlots # enhancements to basic plots
#using CairoMakie # different plotting system
#using AlgebraOfGraphics # even more plotting - similar to ggplot2, uses Makie
using CSV #to read CSV data
using DataFrames #for data wrangling 
#using DataFramesMeta #for data wrangling similar to dplyr functionality
# using ParetoSmooth #for LOO cross-validation computations
# using DifferentialEquations
#using LinearAlgebra: I #identity matrix, for use in models
#using Profile #for profiling code to figure out bottlenecks
#using StatProfilerHTML #nice formatting of profile output
#using RData #to load files in Rdata format

## ---- setup --------
# set random seed
Random.seed!(1234)



## ---- loaddata --------
datafile = "simdat.csv"
#dataFile = "../posts/2022-02-22-longitudinal-multilevel-bayes-1/simdat.Rds"
# load data into a DataFrame object - uses DataFrames package
dat = CSV.read(datafile, DataFrame)


## ---- inspectdata --------
describe(dat)
# set values for priors mu_a, sig_a, mu_b,...
pv = [1.0, 1.0, 1.0, 1.0, 2.5, 1.0, 0.3, 1.0]
outcome = dat.outcome
time = dat.time
Nind = length(unique(dat.id))
Nobs = length(dat.outcome)
id = dat.id # vector keeping track of which data point belongs to which individual

## ---- define-turing-model ----
######################################
## Turing code for model
## pr contain values for prior distributions
######################################
@model function turingmod(outcome, time, pv, Nind, Nobs, id)
    
    ##########
    # priors
    ##########
    # residual population variation
    sig ~ Exponential(1)
    
    # individual variation of each model parameter
    a0 ~ filldist(Normal( pv[1] , pv[2] ),Nind);
    b0 ~ filldist(Normal( pv[3] , pv[4] ),Nind);
    g0 ~ filldist(Normal( pv[5] , pv[6] ),Nind);
    e0 ~ filldist(Normal( pv[7] , pv[8] ),Nind);

    #set up empty variables to be filled below
    # virus_pred = zeros(Nobs) 

    # compute main model parameters
    # note that dos_adj contains a dose for each observation, 
    # but it's repeated after the first Nind entries, so we don't need to loop over all obersvations
        alph = exp.(30*a0)
        bet = exp.(b0)
        gamm = exp.(g0)
        et = exp.(e0)
    
    # loop over all observations
    # since paramaters are saved in vectors of length corresponding to number of individuals
    # we need to index with that extra id[i] notation
    for i in 1:Nobs
      # deterministic model
      n1 = 2*alph[id[i]]
      d1 = exp( -bet[id[i]] * (time[i] - gamm[id[i]]) );
      d2 = exp(   et[id[i]] * (time[i] - gamm[id[i]]) ); 
      virus_pred = log( n1   / ( d1 + d2) );
      outcome[i] ~ Normal(virus_pred, sig)
    end

 end #end turing model block


 turingmod1 = turingmod(outcome, time, pv, Nind, Nobs, id)



## ---- diagnose_mods ----
function modeldiagnostics(postsamp,priorsamp, parnames)
    
    # pull out parameters we want to explore
    # Not sure why this weird syntax works, got it from here:
    # https://discourse.julialang.org/t/debugging-and-visualizing-high-dimensional-models-turing-mcmcchains/41872
    subsamp = postsamp[:,parnames,:]
    subsamppr = priorsamp[:,parnames,:]
    
    # model performance
    modperform = summarystats(subsamp)
    
    # parameter distributions
    pardist = quantile(subsamp)
    
    # trace, density, correlation, mean trace, autocorrelation plots
    tp = StatsPlots.plot(subsamp, seriestype = :traceplot)
    dp = StatsPlots.plot(subsamp, seriestype = :mixeddensity)
    cp = StatsPlots.plot(subsamp, seriestype = :corner)
    mp = StatsPlots.plot(subsamp, seriestype = :meanplot)
    acp = StatsPlots.plot(subsamp, seriestype = :autocorplot)

    # comparing prior and posterior distribution for parameters 
    # convert MCMCChains object returned from Turing into a data frame
    dfpost = DataFrame(subsamp)
    dfprior = DataFrame(subsamppr)

    # select and transform the 2 data frames into a single one
    # uses the Chain and DataFramesMeta packages
    # note that Turing returns objects of type Chains and there is the
    # completely unrelated Chain package for doing R pipe-style coding
    dfpo2 = @chain dfpost begin
        @select $(parnames) # keep only model parameters
        @transform(:type = "post") #add column indicating it's posterior
    end

    dfpr2 = @chain dfprior begin
        @select $(parnames) 
        @transform(:type = "prior")
    end

    df = [dfpr2;dfpo2] #combine data frames

    # turn into long format for plotting
    dflong = DataFrames.stack(df, parnames, [:type], variable_name=:parameter, value_name=:value)

    # prior and posterior density distribution
    # facetted by paremeter
    ppplot = draw( data(dflong) * 
            mapping(:value, color=:type, layout = :parameter) *
                       AlgebraOfGraphics.density() 
                   )

    # put all results in a dictionary (similar to an R list)
    moddiag = Dict("modperform" => modperform, 
                   "pardist" => pardist,
                   "traceplot" => tp,
                   "densityplot" => dp,
                   "corrplot" => cp,
                   "meantraceplot" => mp,
                   "autocorrplot" => acp,
                   "ppplot" => ppplot
                    )

    return moddiag
end #end function that does diagnostics




## ---- fitconditions_m1 ----
# Specify settings for sampler
iterations = 200
warmup = 100
adapt_delta = 0.8
max_td = 5
chains = 2


# Do the sampling/fitting 
res1 = sample(turingmod1, NUTS(warmup, adapt_delta), MCMCThreads(), iterations, chains; max_td)
#res1 = sample(mod1, NUTS(warmup, adapt_delta), MCMCThreads(), iterations, chains, init_params = init_vals; max_td)
# to get prior distributions
res1pr =  sample(turingmod1, Prior(), MCMCThreads(), iterations, chains;)


## ---- prep_explore_m1 ----
postsamp = res1 #consistent naming
priorsamp = res1pr
parnames = names(res1, :parameters) #parameters to show in outputs
# need to do further filtering to get the ID specific parameters removed

## ---- savefits ----
write("m1res_post.jls", postsamp)
write("m1res_prior.jls", priorsamp)

##uncomment to load saved chains
postsamp = read("m1res_post.jls", Chains)
priorsamp = read("m1res_prior.jls", Chains)






## ---- diagnose_m1 ----
## get summary results and diagnostics
summarystats(postsamp)
# seems equivalent
# summarize(res1)
# MCMCChains.summarize(res1; sections = :parameters)

## ---- summarize_m1 ----
# parameter distributions
quantile(postsamp)

## ---- plots_m1 ----
## various plots
# trace, density, correlation, mean trace, autocorrelation
tp = StatsPlots.plot(postsamp, seriestype = :traceplot)
dp = StatsPlots.plot(postsamp, seriestype = :mixeddensity)
cp = StatsPlots.plot(postsamp, seriestype = :corner)
mp = StatsPlots.plot(postsamp, seriestype = :meanplot)
acp = StatsPlots.plot(postsamp, seriestype = :autocorplot)

## ---- priorplots_m1 ----
# comparing prior and posterior distribution for parameters 

# convert MCMCChains object returned from Turing into a data frame
dfpost = DataFrame(postsamp)
dfprior = DataFrame([priorsamp])
parnames = names(res1, :parameters)

# select and transform the 2 data frames into a single one
# uses the Chain and DataFramesMeta packages
# note that Turing returns objects of type Chains and there is the
# completely unrelated Chain package for doing R pipe-style coding
dfpo2 = @chain dfpost begin
    @select $(parnames) # keep only model parameters
    @transform(:type = "post") #add column indicating it's posterior
end

dfpr2 = @chain dfprior begin
    @select $(parnames) 
    @transform(:type = "prior")
end

df = [dfpr2;dfpo2] #combine data frames

# turn into long format for plotting
dflong = DataFrames.stack(df, parnames, [:type], variable_name=:parameter, value_name=:value)

# prior and posterior density distribution
# facetted by paremeter
# uses the AlgebraOfGraphics package 
# (which is meant to be similar to ggplot)
plt = draw(
        data(dflong) * 
        mapping(:value, color=:type, layout = :parameter) *
                     AlgebraOfGraphics.density() 
           )
save(joinpath(figDir,"ppplot_m1.png"), plt)

## ---- preds_m1 ----
## predictive checks
## need to set outcome/Y-values to missing so they will be predicted
Y_new = Vector{Missing}(missing, length(Y)) # vector of `missing`
model_pred = mod1turing(C, Y_new)
pred_post = predict(model_pred, res1)  # posterior
pred_prior = predict(model_pred, res1pr)

# distributions for prior/posterior predicted outcomes
postpred = quantile(pred_post)
priorpred = quantile(pred_prior)

## ---- make_predplot1_m1 ----
# plot showing prior and posterior predictions
f = Figure()
Axis(f[1, 1], xlabel = "Concentration", ylabel = "Reduction", title = "Prior predictive check")
CairoMakie.scatter!(dat2.C, dat2.Y)
# 4th column is 50% quartile
CairoMakie.lines!(dat2.C, priorpred[:, 4])
# column 2 and 6 are 2.5 and 97.5 quartiles
CairoMakie.band!(
    dat2.C,
    priorpred[:, 2],
    priorpred[:, 6],
    color = (:blue, 0.2),
)
Axis( f[1, 2], xlabel = "Concentration", ylabel = "Reduction", title = "Posterior predictive check")
CairoMakie.scatter!(dat2.C, dat2.Y)
# 4th column is 50% quartile
CairoMakie.lines!(dat2.C, postpred[:, 4])
# column 2 and 6 are 2.5 and 97.5 quartiles
CairoMakie.band!(
    dat2.C,
    postpred[:, 2],
    postpred[:, 6],
    color = (:blue, 0.2),
)
save(joinpath(figDir,"pppred_m1.png"), f)


