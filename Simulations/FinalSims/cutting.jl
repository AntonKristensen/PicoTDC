using CSV
using DataFrames


function cut(data, cutparameterfile = "..cutparams.csv")
    cutread = CSV.read(cutparameterfile, DataFrame; header=false, delim=",", ignorerepeated=false)
    cutparams = collect(cutread[1,:])

    cutfunc(E, p) = p[1] .+ p[2] * exp.(- (E .+ p[3]) ./ p[4])

    cutlee = 2 # The cut isn't perfect, so add a bit of extra leeway
    lowert = 0.5 # Set a lower threshold for detection. Usefull for cutting out crosstalk also

    incidents = data[:,1]
    firsts = data[:,2]
    seconds = data[:,3]
    #fronts = data[:,4]
    #backs = data[:,5]

    cutindices = firsts .< cut(incidents, cutparams) .+ cutlee .&& firsts .> lowert .&& seconds .> lowert .&&  seconds .< cut(incidents, cutparams) .+ cutlee
    
    return data[cutindices, :]

end

