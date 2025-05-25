

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


include("functions.jl")

@time collected = collector("output/back1.phsp")

print(collected)

#incidentframe = detectorlooping(minimumenergy=1, addedtime = 1e-6)
#matchwriter(incidentframe)







