

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


include("functions.jl")


incidentframe = detectorlooping(minimumenergy=1, addedtime = 1e-6)
matchwriter(incidentframe)







