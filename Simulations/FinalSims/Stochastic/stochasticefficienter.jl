
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Glob

include("../cutting.jl")


files = glob("output/matches*")
for filepath in files
    noncutdata = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)
    cleaneddata = noncutdata[noncutdata[:,1] .!= "Incident energy", :] # For some reason the matchwriter doesn't see that the file isn't empty, so it writes the header again in the middle of a lot of data :/
    println(noncutdata)
    data = cut(cleaneddata)

    ###############
    # Collecting results from all individual pairs into one big list
    incidents = data[:,1]
    firsts = data[:,2]
    seconds = data[:,3]
    fronts = data[:,4]
    backs = data[:,5]
    frontevents = data[:,6]
    backevents = data[:,7]


    cincidents = incidents[frontevents .== backevents]
    fincidents = incidents[frontevents .!= backevents]


    println("Correct: ", length(incidents[frontevents .== backevents]), ", False: ", length(incidents[frontevents .!= backevents]))


    geometry = CSV.read("geometry.txt", DataFrame; delim=",", ignorerepeated=true, ignoreemptyrows=true)

    # Making it so that the neutron beam is wide enough to cover all the front scintillators. 
    frontmaxx = string(maximum(abs.(geometry[geometry[:,end], 1]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
    frontmaxy = string(maximum(abs.(geometry[geometry[:,end], 2]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
    frontmaxsize = string(maximum(abs.(geometry[geometry[:,end], 4]./2)))


    correct = length(incidents[frontevents .== backevents])
    fake =  length(incidents[frontevents .!= backevents])

    println("Correct: ", correct, ", Fake: ", fake)
end
