

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


#include("../functions.jl")
include("../testfunctions.jl")

incidentframe = detectorlooping(timeuncertainty=5e-11, minimumenergy=1, addedtime = 0.02e-9) # DPCT shoots roughly 50 protons/nanosecond
matchwriter(incidentframe)







