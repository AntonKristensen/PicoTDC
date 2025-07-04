

using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Glob



include("../testfunctions.jl")

#times = [100 333 1000 3333 10000 33333 100000 333333 1000000]
times = [10000 33333 100000 333333 100000]
#times = collect(10000:5000:50000)

for t in times
    println(t)
    @time incidentframe = detectorlooping(minimumenergy=1, addedtime = t * 1.0e-12, threshold=0.25)
    matchwriter(incidentframe, file = "output/matches" * string(t) * ".csv")
end





