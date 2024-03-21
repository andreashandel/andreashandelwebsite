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
# load/install all needed packages
using CSV #to read CSV data
using DataFrames #for data wrangling 
using DataFramesMeta #for data wrangling similar to dplyr functionality
using Random #for random numbers
using Distributions #for probabilistic distributions
using Turing #Bayesian fitting routines
using Chain #for chaining commands like the R pipe
using MCMCChains #for visualization of turing output
using Plots # for basic plots
using StatsPlots # enhancements to basic plots
using CairoMakie # different plotting system
using AlgebraOfGraphics # even more plotting - similar to ggplot2, uses Makie
using ParetoSmooth #for LOO cross-validation computations
using DifferentialEquations
#using LinearAlgebra: I #identity matrix, for use in models
#using Profile #for profiling code to figure out bottlenecks
#using StatProfilerHTML #nice formatting of profile output
#using RData #to load files in Rdata format



## ---- setpaths --------
######################################
# Setting paths and file names
# using the DrWatson functions is similar to here() in R
######################################
# set random seed
Random.seed!(1234)



## ---- loaddata --------
datafile = "simdat.csv"
#dataFile = "../posts/2022-02-22-longitudinal-multilevel-bayes-1/simdat.Rds"
# load data into a DataFrame object
# this uses the RData package
# note that the file needs to end in .rds for some reason .Rds doesn't work, I had to rename
# https://github.com/JuliaData/RData.jl/issues/91
dat = CSV.read(datafile, DataFrame)


## ---- inspectdata --------
describe(dat)




## ---- defineode --------
# original model
function odemod1!(dy, y, p, t)
    alph, bet, gamm, et = p
    dy[1] = - bet*y[1]*y[3] 
    dy[2] = bet*y[1]*y[3] - gamm*y[2]
    dy[3] = alph*y[2]  - et*y[3]
end

# same model as above, but simulated in log space
function odemod2!(dy, y, p, t)
    alph, bet, gamm, et = p
    dy[1] = - bet * exp(y[3]);
    dy[2] = bet * exp(y[1] + y[3] - y[2]) - gamm;
    dy[3] = alph * exp(y[2] - y[3]) - et;
end


## ---- runode --------
y0 = [log(1e8), 0, 5]
tspan = (0.0, 10.0)
p = [3.0, -1.0, 1.0, 1.0]
# run model 1
odeprob1 = ODEProblem(odemod1!, y0, tspan, p)
odesol1 = solve(odeprob1)
Plots.plot(odesol1, layout = (3,1))
# run model 2
odeprob2 = ODEProblem(odemod2!, log(y0), tspan, p)
odesol2 = solve(odeprob2)
Plots.plot(odesol2, layout = (3,1))




## ---- define-turing-model ----
######################################
## Turing code for model
## pr contain values for prior distributions
######################################
@model function turingmod(pr, outcome, time, )
    
    ##########
    # priors
    ##########
    # residual population variation
    sig ~ Exponential(1)
    
    # hyper-priors to allow for adaptive pooling among individuals 
    mu_a ~ Normal(pr["mu_a_mu"], pr["mu_a_sd"])
    mu_b ~ Normal(pr["mu_b_mu"], pr["mu_b_sd"])
    sig_a ~ Exponential(1)
    sig_b ~ Exponential(1)
    
    # individual variation of each model parameter
    a0 ~ filldist(Normal( mu_a , sig_a ),Nind);
    b0 ~ filldist(Normal( mu_b , sig_b ),Nind);
  
    # average dose-dependence of each ODE model parameter
    a1 ~ Normal( pr["a1_mu"] , pr["a1_sd"]); 
    b1 ~ Normal( pr["b1_mu"] , pr["b1_sd"]);


    # compute main model parameters
    # note that dos_adj contains a dose for each observation, 
    # but it's repeated after the first Nind entries, so we don't need to loop over all obersvations
    for  i in 1:Nind 
        alpha[i] = a0[i] + a1 * dose_adj[i]
        beta[i] = b0[i] + b1 * dose_adj[i]
    end
    
    # loop over all observations
    # since paramaters are saved in vectors of length corresponding to number of individuals
    # we need to index with that extra id[i] notation
    for i in 1:Nobs
      # deterministic model
      virus_pred[i] = exp(alpha[id[i]]) * log(time[i]) - exp(beta[id[i]]) * time[i] 
      # likelihood
      Y[i] ~ Normal(virus_pred[i], sig)
    end


 end #end model block

bayesmod = turingmod(pr,outcome,time,Nobs,Nind)
