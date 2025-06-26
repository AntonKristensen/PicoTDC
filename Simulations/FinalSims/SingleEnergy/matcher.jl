

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


#include("../functions.jl")
include("../testfunctions.jl")

incidentframe = detectorlooping(minimumenergy=1, addedtime = 1e-6)
matchwriter(incidentframe)

thirtyuncertainframe = detectorlooping(timeuncertainty=3e-11, minimumenergy=1, addedtime = 1e-6)
matchwriter(thirtyuncertainframe, file="output/thirtyuncertainmatches.csv")

fiddiuncertainframe = detectorlooping(timeuncertainty=5e-11, minimumenergy=1, addedtime = 1e-6)
matchwriter(fiddiuncertainframe, file="output/fiddiuncertainmatches.csv")

hunniuncertainframe = detectorlooping(timeuncertainty=1e-10, minimumenergy=1, addedtime = 1e-6)
matchwriter(hunniuncertainframe, file="output/hunniuncertainmatches.csv")







