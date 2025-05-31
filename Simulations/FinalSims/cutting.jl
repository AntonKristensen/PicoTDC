using CSV
using DataFrames


function cut(data, cutparameterfile = "../cutparams.csv")
    
    cutread = CSV.read(cutparameterfile, DataFrame; header=false, delim=",", ignorerepeated=false)
    cutparams = collect(cutread[1,:])
	
    cutlee = 1.5 # The cut isn't perfect, so add a bit of extra leeway
    lowert = 0.5 # Set a lower threshold for detection. Usefull for cutting out crosstalk also
    minimumenergy = 25 # Sets an energy below which the cut will not be applied

    #cutfunc(E, p) = p[1] .+ p[2] * exp.(- (E .+ p[3]) ./ p[4]) # Old function, exponentially decaying. Not good at describing the behaviour.
    cutfunc(E, p) = p[1] .+ p[2] ./ (E .+ p[3])

    incidents = data[:,1]
    firsts = data[:,2]
    seconds = data[:,3]
    #fronts = data[:,4]
    #backs = data[:,5]

    cutindices = firsts .< cutfunc(incidents, cutparams) .+ cutlee .&& firsts .> lowert .&& seconds .> lowert .&&  seconds .< cutfunc(incidents, cutparams) .+ cutlee
    
    return data[cutindices .|| incidents .< minimumenergy, :]

end

