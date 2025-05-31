

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


stochastictime = parse(Float64, ARGS[1]) * 1e-12 # The added time is put in to the julia program as argument, with units of picoseconds

println(stochastictime, " seconds")

include("../testfunctions.jl")

println("Tråde: ", Threads.nthreads())

@time incidentframe = detectorlooping(minimumenergy=1, addedtime = stochastictime, threshold=0.25)
matchwriter(incidentframe)







