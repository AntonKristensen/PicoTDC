
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
function tracereader(file, channel=0)
    lines = readlines(file)

    traces = []
    for i in 1:length(lines)
        if length(lines[i]) > 20
            if lines[i][1:20] == " <trace channel=\""*string(channel)*"\">"
                append!(traces, [parse.(Int, split(lines[i][21:end]))])
            end
        end

    end
    return traces
end

channel0 = []
channel1 = []
files = glob("*","data/")
for filepath in files
    append!(channel0, [tracereader(filepath, 0)])
    append!(channel1, [tracereader(filepath, 1)])
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

function areafunction(trace, peakindex, threshold=10)
    baseline = mean(trace[1:peakindex]) # Robustly guessing the baseline
    walkindex = peakindex # Start at the peak, then walk back until the pulse starts
    while trace[walkindex] > baseline + threshold
        walkindex -= 1
    end
    area = trace[walkindex] - baseline # The first point starts under the threshold, but probably doesn't contribute a lot to the area
    walkindex += 1
    while trace[walkindex] > baseline + threshold
        area += trace[walkindex] - baseline
        walkindex += 1
    end
    return area
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
    areas = []
    if !isempty(indices) 
        indices, proms = peakproms(indices, trace, min=minprom)
        if !isempty(indices)
            indices, widths = peakwidths(indices, trace, proms)            
            heights = trace[indices] # The last on here is the height of the peaks after filtering
            
            for j in 1:length(indices) # Finding peak position and heights by fitting 2nd order polynomial on top 3 points, and areas by my algorithm
                fitindex, fitheight = heightfunction(trace, indices[j])
                append!(fitindices, fitindex)
                append!(fitheights, fitheight)
                append!(areas, areafunction(trace, indices[j]))
            end

        end
    end
    intervals = [] # Finding intervals between peaks
    if length(indices) > 1


        append!(intervals, fitindices[2:end] - fitindices[1:end-1])
    end
    return indices, proms, widths, heights, maximum(trace), intervals, fitindices, fitheights, areas # The maximum is of the entire trace, so it is a single value, not a vector
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




###########################################################

#coincplot = plot()
println(files)
println(length.(channel0))
filenumber = 1
for n in 1:length(channel1[filenumber])
    pulses = pulseextractor(channel1[filenumber][n], minprom = 10, minheight=2200)
    if pulses[5] > 2250 # If the maximum value of the trace is high enough, do analysis
        println(pulses)
        coincplot = plot()
        plot!(channel0[filenumber][n], color=:blue, alpha=0.5, label="First")
        plot!(channel1[filenumber][n], color=:red, alpha=0.5, label="Second")
        title!(files[filenumber])
        display(coincplot)
    end
end
#display(coincplot)
