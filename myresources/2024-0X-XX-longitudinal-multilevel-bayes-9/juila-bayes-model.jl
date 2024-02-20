########################################################
# Model fitting with Turing/Julia 
# Fitting longitudinal data with a Bayesian hierarchical model
# Last updated 2024-02-16 (or later) by Andreas Handel
########################################################



######################################
# set up project with DrWatson package
######################################

######################################
# Initializing project
# only needs to happen once
######################################
# set up folder structure for project
# could use the defaults, but we want our own structure
# this loosely follows the MRG msproject template
# folders = [
#         "data" => ["source", "derived"],
#         "deliv", => ["figure", "table", "report"],
#         "docs",
#         "model",
#         "script"
#         ]
#DrWatson.initialize_project("ex2-turing-julia"; template = folders, git = false)


## ---- setup --------
######################################
# Loading project
# should happen each time
######################################
# load DrWatson
# other packages are loaded once project is activated
using DrWatson
# those 2 seem equivalent
DrWatson.quickactivate(@__DIR__)
#DrWatson.quickactivate("ex2-julia")


## ---- loadpackages --------
# now that project is activated, load all needed packages
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
#using LinearAlgebra: I #identity matrix, for use in models
#using Profile #for profiling code to figure out bottlenecks
#using StatProfilerHTML #nice formatting of profile output

