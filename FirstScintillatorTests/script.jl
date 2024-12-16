
using Plots
using CSV
using DataFrames
using Glob


################################### Function definitions
function fwhm(x,y)
    value, index = findmax(y)

    leftindex = 1
    while y[leftindex] < value/2
        leftindex += 1
    end

    rightindex = length(x)
    while y[rightindex] < value/2
        rightindex -= 1
    end

    return x[rightindex] - x[leftindex]
end

function risetime(x,y)
    value, index = findmax(y)

    start = 1
    while y[start] < value/10
        start += 1
    end
    stop = start
    while y[stop] < value*9/10
        stop += 1
    end
    return x[stop] - x[start]
end

function pulseplot(file)
    data = CSV.read(file, DataFrame, header=12)
    times = data[:,1]  *1e7
    volts = data[:,2]
    plottinate = plot(times,volts)
    display(plottinate)
    xlabel!("Time")
end
##################################################

nor = glob("0OhmTests/*csv")

bigr = glob("100OhmTests/*csv")

noparams = zeros(3,length(nor))
for n in 1:length(nor)
    data = CSV.read(nor[n], DataFrame, header=12)
    times = data[:,1]*1e7
    volts = data[:,2]
    noparams[1,n] = fwhm(times,volts) #fwhm
    noparams[2,n] = sum(volts)        #area (charge)
    noparams[3,n] = risetime(times,volts)                 #rise time
end

bigparams = zeros(3,length(bigr))
for n in 1:length(bigr)
    data = CSV.read(bigr[n], DataFrame, header=12)
    times = data[:,1]*1e7
    volts = data[:,2]
    bigparams[1,n] = fwhm(times,volts) #fwhm
    bigparams[2,n] = sum(volts)        #area (charge)
    bigparams[3,n] = risetime(times,volts)              #rise time
end

fwhmplot = scatter(noparams[2,:],noparams[1,:],seriestype=:scatter,color=:blue,label="0 Ohm")
fwhmplot = scatter!(bigparams[2,:],bigparams[1,:],seriestype=:scatter,color=:red,label="100 Ohm")
title!("FWHM vs area of pulse")
xlabel!("Pulse area")
ylabel!("FWHM (ns)")
display(fwhmplot)

risetimeplot = scatter(noparams[2,:],noparams[3,:],seriestype=:scatter,color=:blue,label="0 Ohm")
risetimeplot = scatter!(bigparams[2,:],bigparams[3,:],seriestype=:scatter,color=:red,label="100 Ohm")
title!("Risetime vs area of pulse")
xlabel!("Pulse area")
ylabel!("Risetime (ns)")
display(risetimeplot)

bandwidthplot = scatter(noparams[2,:],0.35 ./ noparams[3,:],seriestype=:scatter,color=:blue,label="0 Ohm")
bandwidthplot = scatter!(bigparams[2,:],0.35 ./bigparams[3,:],seriestype=:scatter,color=:red,label="100 Ohm")
title!("Bandwidth vs area of pulse")
xlabel!("Pulse area")
ylabel!("Bandwitdh (GHz)")
display(bandwidthplot)
