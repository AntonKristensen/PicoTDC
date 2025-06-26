

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


#include("../functions.jl")
include("../testfunctions.jl")

println("Tråde: ", Threads.nthreads())

@time incidentframe = detectorlooping(minimumenergy=1, addedtime = 1e-6, threshold=0.25)
matchwriter(incidentframe)

@time fiddiincidentframe = detectorlooping(timeuncertainty=5e-11, minimumenergy=1, addedtime = 1e-6, threshold=0.25)
matchwriter(fiddiincidentframe, file="output/fiddiuncertainmatches.csv")





