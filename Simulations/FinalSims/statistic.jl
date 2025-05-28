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
    medi = mode
    lower = medi - quantile(data, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
    upper = quantile(data, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak
    println("Brute: ", medi, ", σ-: ", lower, ", σ+:", upper)

    bound = min(lower,upper) # It works because the right side of the peak is quite nicely gaussian!

    # ML fit, cutting data 2 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound*2) .& (data .< medi + bound*2)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    println("Fit:", params(gaussfit))



    return params(gaussfit)
end

