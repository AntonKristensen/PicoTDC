
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)



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

function statisticing(data)
    # Doing a slight bit of statistics
    medi = median(data) # Getting a decent robust guess for the mean of the peak so I can make a un-bad cut when fitting
    lower = medi - quantile(data, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
    upper = quantile(data, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak

    #print(medi, ", ", lower, ", ", upper, "\n")

    bound = min(lower,upper)

    # ML fit, cutting data 2 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound*2) .& (data .< medi + bound*2)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    #print(medi, " ", bound, " ", params(gaussfit), "\n")
    return params(gaussfit)
end


println("Correct: ", length(incidents[frontevents .== backevents]), ", False: ", length(incidents[frontevents .!= backevents]))


geometry = CSV.read("geometry.txt", DataFrame; delim=",", ignorerepeated=true, ignoreemptyrows=true)

# Making it so that the neutron beam is wide enough to cover all the front scintillators. 
frontmaxx = string(maximum(abs.(geometry[geometry[:,end], 1]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
frontmaxy = string(maximum(abs.(geometry[geometry[:,end], 2]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
frontmaxsize = string(maximum(abs.(geometry[geometry[:,end], 4]./2)))

correct = length(incidents[frontevents .== backevents])
fake =  length(incidents[frontevents .!= backevents])
neutrons = parse(Float64, ARGS[1]) * parse(Float64, ARGS[2])
flux = neutrons / (parse(Float64, frontmaxx) * parse(Float64, frontmaxy))
matchrate = correct / (flux / 10000)
println(matchrate, " matches / (neutrons/cm²)")

writingfile = open("geometryresults.txt", "a")
write(writingfile, string(matchrate) * "\n")
close(writingfile)
