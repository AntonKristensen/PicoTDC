using Plots
using Glob
using DelimitedFiles
using Statistics
using Distributions
using StatsBase
using LsqFit

function nodip(t, p)
    booler = t .> 0-p[4]
    return (p[3] * (exp.(-(t .+p[4])/p[2]) .- exp.(-(t .+p[4])/p[1])) .* booler) .+ p[5]
end

function doublexp(t, p) # p is the parameters, where p[1] is the rise time, p[2] is the decay time and p[3] is the amplitude and p[4] is the time shift
    booler = t .> 0-p[4]
    return p[3] * (exp.(-(t .+p[4])/p[2]) .- exp.(-(t .+p[4])/p[1])) .* booler
end

function doubledipper(t, p)
    return (doublexp(t,p[1:4]) .- doublexp(t,[p[5], p[6], p[7], p[4]])) .+ p[8]
end

filenames = glob("*.txt","SiPM test Anton/")

fig = plot(legend=false,dpi=600,minorgrid = true, minorticks=4,minorgridstyle=:solid)
risetimes = []
fwhms = []
decaytimes = []
for file in filenames
    data = readdlm(file,skipstart=6)

    max = maximum(data[:,2])

    risestart = data[findfirst(x -> x > max /10, data[:,2]),1]
    risestop = data[findfirst(x -> x > max *9 / 10, data[:,2]),1]
    time = risestop-risestart
    push!(risetimes, time)

    decaystop = data[findlast(x -> x > max /10, data[:,2]),1]
    decaystart = data[findlast(x -> x > max *9 / 10, data[:,2]),1]
    time = decaystop-decaystart
    push!(decaytimes, time)

    fwhmstart = data[findfirst(x -> x > max /2, data[:,2]),1]
    fwhmstop = data[findlast(x -> x > max /2, data[:,2]),1]
    push!(fwhms,fwhmstop-fwhmstart)

    plot!(data[:,1],data[:,2],alpha = 0.8)

    #=
    fitguess = [time, 2 * time, 5*max, 1.0e-9, 3*time, 4 * time, 5*max, 1.0e-5]
    fitt = curve_fit(doubledipper, data[:,1], data[:,2], fitguess)
    plot!(data[:,1], doubledipper(data[:,1], fitt.param), linestyle=:dash, linecolor=:black, alpha=0.5)
    #plot!(data[:,1], doubledipper(data[:,1], fitguess), linestyle=:dash, linecolor=:black, alpha=0.5)
    println(fitguess)
    #println(fitt.param)
    display(fig)
    =#


end

println("Risetime = ", round(mean(risetimes),sigdigits=3),"ns ± ", round(std(risetimes),sigdigits=3) / sqrt(length(risetimes)),"ns. σ = ",round(std(risetimes),sigdigits=3) ,"\n")
println("Decaytime = ", round(mean(decaytimes),sigdigits=3),"ns ± ", round(std(decaytimes),sigdigits=3) / sqrt(length(decaytimes)),"ns. σ = ",round(std(decaytimes),sigdigits=3) ,"\n")
println("FWHM = ",round(mean(fwhms),sigdigits=3),"ns ± ", round(std(fwhms) / sqrt(length(fwhms)),sigdigits=3),"ns. σ = ",round(std(fwhms),sigdigits=3) ,"\n")




xlabel!("Time (s)")
ylabel!("Voltage (V)")
title!("Fast output pulses")
xlims!(-0.00000001,0.00000004)
display(fig)
savefig(fig, "pulsesplot.svg")

# Calculate the error in the estimate of pulse energy by ToT method
function Eerror(thresholdratio, dthreshold, dtot, tau, dtau) # Thresholdratio is E/V, where E is the peak voltage of a pulse and V is the voltage threshold for ToT. dthreshold is dV/V so the relative uncertainty of the threshold
    energyratio = 1 ./ thresholdratio
    tot = tau .* log.(energyratio)

    return sqrt.(dthreshold.^2 .+ log.(energyratio).^2 .* ((dtot ./ tot).^2 .+ (dtau ./tau).^2))
end


ratios = collect(0.0:0.001:0.2)
fig2 = plot(ratios,100* Eerror(ratios, 0.01, 0.1, mean(decaytimes)*1e9, std(decaytimes)*1e9), legend=false, linewidth=3)
title!("Pulse Energy uncertainty")
xlabel!("Threshold / pulse max")
ylabel!("Energy uncertainty (%)")
savefig("EnergyError.svg")
display(fig2)

