
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Glob
using Peaks
using StatsBase

# FUNCTIONS
######################

# Function for loading a file from caenscope and returning the traces in a list
function tracereader(file)
    lines = readlines(file)

    traces = []
    for i in 1:length(lines)
        if length(lines[i]) > 20
            if lines[i][1:20] == " <trace channel=\"0\">"
                append!(traces, [parse.(Int, split(lines[i][21:end]))])
            end
        end

    end
    return traces
end

data = []
files = glob("*","data/")
for filepath in files
    append!(data, [tracereader(filepath)])
end

function heightfunction(trace, peakindex) # Finding the height (and position) of a 2nd degree polynomial from 3 points
    # Building Vandermonde matrix
    x1 = peakindex -1
    x2 = peakindex
    x3 = peakindex +1
    V = [1 x1 x1^2 ; 1 x2 x2^2 ; 1 x3 x3^2]
    y = [trace[x1], trace[x2], trace[x3]]
    coefficients = inv(V) * y
    topx = - coefficients[2] / (2 * coefficients[3]) # x_top  = - b / 2a for ax² bx + c = y
    topvalue = coefficients[1] + coefficients[2] * topx + coefficients[3] * topx^2
    return topx, topvalue
end

function pulseextractor(trace; minheight = 0, minprom=0, window=1) # Set a minimum prominence to filter peaks that are just noise
    indices, heights = findmaxima(trace, window)
    indices = indices[heights .> minheight] # Cutting out peaks that are too low
    heights = heights[heights .> minheight]
    widths = []
    heights = []
    proms = []
    fitindices = []
    fitheights = []
    if !isempty(indices) 
        indices, proms = peakproms(indices, trace, min=minprom)
        if !isempty(indices)
            indices, widths = peakwidths(indices, trace, proms)            
            heights = trace[indices] # The last on here is the height of the peaks after filtering
            
            for j in 1:length(indices) # Finding peak position and heights by fitting 2nd order polynomial on top 3 points
                fitindex, fitheight = heightfunction(trace, indices[j])
                append!(fitindices, fitindex)
                append!(fitheights, fitheight)
            end

        end
    end
    intervals = [] # Finding intervals between peaks
    if length(indices) > 1


        append!(intervals, fitindices[2:end] - fitindices[1:end-1])
    end
    return indices, proms, widths, heights, maximum(trace), intervals, fitindices, fitheights # The maximum is of the entire trace, so it is a single value, not a vector
end



# Plot all the traces on top of each other, with option for plotting only the first N traces
function traceplotter(traces, number=length(traces))
    fig1 = plot()
    for n in 1:minimum([length(traces),number])
        plot!(traces[n],alpha=0.2)
    end
    title!("Pulses")
    xlabel!("Time (ADC units ~ 4ns)")
    ylabel!("Voltage (ADC units)")
    display(fig1)
    savefig("plots/traces.png")
end




# DOING STUFF
########################

function dataanalyzer(datasetnumber; plotting=false)
    traces = data[datasetnumber]

    if plotting
        testtrace = traces[2901]
        peakindices, foo = findmaxima(testtrace)
        peakplot = plotpeaks(testtrace; peaks=peakindices, prominences=true, widths=true)
        savefig("plots/traceplot.png")
        display(peakplot)
    end

    indices = []
    prominences = []
    widths = []
    heights = []
    maxs = []
    intervals = []
    fitindices = []
    fitheights = []
    for n in 1:length(traces)
        pulses = pulseextractor(traces[n], minheight = 0, minprom=10, window=1)
        if !isempty(pulses) # Don't do anything if there aren't pulses in the trace (happens if minimum prominence is more than the trigger level)
            append!(indices, pulses[1])
            append!(prominences, pulses[2])
            append!(widths, pulses[3])
            append!(heights, pulses[4])
            append!(maxs, pulses[5])
            append!(intervals, pulses[6])
            append!(fitindices, pulses[7])
            append!(fitheights, pulses[8])
        end
    end


    # PLOTTING!
    #####################
    if plotting
        promplot = histogram(prominences)
        title!(files[datasetnumber][6:end-4])
        xlabel!("Pulse prominences")
        ylabel!("Counts")
        savefig("plots/prominencehist.png")
        display(promplot)
        sleep(1)
        heightplot = histogram(heights)
        title!(files[datasetnumber][6:end-4])
        xlabel!("Pulse heights")
        ylabel!("Counts")
        savefig("plots/heighthist.png")
        display(heightplot)
        sleep(1)
        widthplot = histogram(widths[widths .< 10])
        title!(files[datasetnumber][6:end-4])
        xlabel!("Pulse Widths")
        ylabel!("Counts")
        savefig("plots/widthhist.png")
        display(widthplot)
        sleep(1)
        maxplot = histogram(maxs, bins=500)
        title!(files[datasetnumber][6:end-4])
        xlabel!("Trace maxima")
        ylabel!("Counts")
        savefig("plots/maxhist.png")
        display(maxplot)
        sleep(1)
        intervalplot = histogram(intervals, bins=1000)
        title!(files[datasetnumber][6:end-4])
        xlabel!("Peak intervals")
        ylabel!("Counts")
        savefig("plots/peakintervals.png")
        display(intervalplot)
        sleep(1)

        doublehistogram = histogram2d(widths[widths .< 10], heights[widths .< 10], bins=(500,500))
        title!(files[datasetnumber][6:end-4])
        xlabel!("Pulse width")
        ylabel!("Pulse height")
        savefig("plots/doublehistogram.png")
        display(doublehistogram)
    end


    ######### Calculating on histogram, to get how much time passes between each cyclotron pulse
    edg = 0:maximum(intervals)/1000:maximum(intervals)
    midpoints = collect(0:maximum(intervals)/1000:maximum(intervals) * 999 / 1000) .+ maximum(intervals)/1000
    h = fit(Histogram, intervals, edg)
    bincounts = h.weights

    intervalindices, intervalheights = findmaxima(h.weights, 5)
    
    if plotting
        pulseintervalplot = plot(midpoints,bincounts, color=:blue)
        scatter!(midpoints[intervalindices], bincounts[intervalindices], color=:red)
        savefig("plots/pulseintervals.png")
        display(pulseintervalplot)
    end

    positions = []
    for i in 1:length(intervalindices)
        position, value = heightfunction(bincounts, intervalindices[i])
        append!(positions, position)
    end
    timeintervals = 4 * (positions[2:end] - positions[1:end-1])* maximum(intervals)/1000
    filteredtimeintervals = timeintervals[timeintervals .> 10 .&& timeintervals .< 20]

    if plotting
        timehist = histogram(filteredtimeintervals, bins=50)
        savefig("plots/timehist.png")
        display(timehist)
    end

    #println("Cyclotron interval is: ", mean(filteredtimeintervals), " ± ", std(filteredtimeintervals))

    return mean(filteredtimeintervals), std(filteredtimeintervals)/sqrt(length(filteredtimeintervals))
end

results = zeros(length(files))
errors = zeros(length(files))
datanames = Vector{String}(undef, length(files))
Threads.@threads for i in 1:length(files) # Analyzing all the data sets multithreadedly
    r,e = dataanalyzer(i)
    results[i] = r
    errors[i] = e
    datanames[i] = files[i][6:end-4]
end


resultplot = scatter(datanames, results, yerr = errors, msc = 1, label = "Data", xrotation = 45)
weightedmean = mean(results, weights(1 ./ (errors).^2))
weighteduncertainty = sqrt(1 / sum(1 ./ (errors).^2))
println("Cyclotron frequency: ", weightedmean, " ± ", weighteduncertainty, "ns")
hline!([weightedmean ], color=:red, label="Mean = " * string(round(weightedmean, sigdigits=3)) * "ns")
hspan!([weightedmean - weighteduncertainty, weightedmean + weighteduncertainty], color=:red, alpha=0.3, label="± " * string(round(weighteduncertainty, sigdigits=3)) * "ns")
title!("Cyclotron times from first experiment")
xlabel!("Measurement number")
ylabel!("Cyclotron pulse time (ns)")
savefig("plots/resultplot.png")
display(resultplot)

chi = sum(((results .- weightedmean)./errors).^2)
print("Chi = ", chi, ", DoF = ", length(results) -1)

dataanalyzer(6, plotting=true)

