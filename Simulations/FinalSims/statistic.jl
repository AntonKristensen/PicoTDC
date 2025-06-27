using Plots
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions

function statisticing(data)
    # Doing a slight bit of statistics
    
    # Putting it into a histogram to find the max value
    h = fit(Histogram, data, nbins=100) 
    edges = collect(h.edges[1])
    maxindex = argmax(h.weights)
    mode = (edges[maxindex]+ edges[maxindex+1])/2 # Turns out to be decently robust

    #medi = median(data[data .> 80 .&& data .< 120]) # Getting a decent robust guess for the mean of the peak so I can make a un-bad cut when fitting
    medi = maximum([mode, median(h.weights)])  
println(medi)
    #lower = medi - quantile(data, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
    #upper = quantile(data, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak
    #println("Brute: ", medi, ", σ-: ", lower, ", σ+:", upper)

    #bound = min(lower,upper) # It works because the right side of the peak is quite nicely gaussian!

    bound = medi * 0.15
    # ML fit, cutting data 2 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound) .& (data .< medi + bound)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    println("Fit:", params(gaussfit), ", Data between: ", medi-bound, " to ", medi+bound)



    return params(gaussfit)
end

function resampler(data)
    newindices = rand(1:length(data), length(data))
    return data[newindices]
end

function bootstatisticing(data, bootnumber=10000) # Function for finding mean and spread, but gives uncertainty estimates from bootstrapping
    # Doing a slight bit of statistics
    
    # Putting it into a histogram to find the max value
    h = fit(Histogram, data, nbins=100) 
    edges = collect(h.edges[1])
    maxindex = argmax(h.weights)
    mode = (edges[maxindex]+ edges[maxindex+1])/2 # Turns out to be decently robust
    println("mode: ", mode)
   
    centguess = maximum([mode, median(h.weights)])
println("Center guess:", centguess)


    boundsfit = fit_mle(Normal, data[data .> mode * 0.90 .&& data .< mode * 1.10]) # Fitting a gaussian on only the right side of the data to get a decent idea of where to cut it
    println("Bounds fit: ", boundsfit)
    center = params(boundsfit)[1]
    bound = params(boundsfit)[2] # For a 1 sigma region

    # ML fit, cutting data 1 sigma below and above the calculated mean
    fitdata = data[(data .> center - bound) .& (data .< center + bound)] # Cutting a roughly 3sigma region around the peak

    means = zeros(bootnumber)
    spreads = zeros(bootnumber)
    println("ping")
    for n in 1:bootnumber
        gaussfit = fit_mle(Normal, resampler(fitdata))
        means[n] = params(gaussfit)[1]
        spreads[n] = params(gaussfit)[2]
    end

    println("Fit:", mean(means), ", ", mean(spreads), ", Data between: ", center-bound, " to ", center+bound)

    # Returns 4 variables which are all 2-element lists: mean&uncertainty, deviation&uncertainty, fitpoints&notfitpoints, selectionmean&bound
    return [median(means), std(means)], [median(spreads), std(spreads)], [length(fitdata), length(data)-length(fitdata)], [center, bound]
end
