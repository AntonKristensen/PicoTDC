
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

fiddidata = CSV.read("output/fiddiuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)
fiddiincidents = fiddidata[:,1]

hunnidata = CSV.read("output/hunniuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)
hunniincidents = hunnidata[:,1]


include("../cutting.jl")
cutdata = cut(data)




###############
# Collecting results from all individual pairs into one big list
incidents = data[:,1]
firsts = data[:,2]
seconds = data[:,3]
fronts = data[:,4]
backs = data[:,5]


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
    print(medi, " ", bound, " ", params(gaussfit), "\n")

    # ML fit, cutting data 1 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound*2) .& (data .< medi + bound*2)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    print(medi, ", ", lower, ", ", upper, " ", params(gaussfit),"\n")

    return params(gaussfit)
end

medi, spread = statisticing(incidents)
fiddimedi, fiddispread = statisticing(fiddiincidents)
hunnimedi, hunnispread = statisticing(hunniincidents)



fig1 = histogram(data[data[:,1] .< medi + spread*2, 1], bins=0:1: medi + spread*2, color=:black, label=:none, alpha=1, size=(500,300), dpi=1000)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/NoCutEnergies.svg")


fig2 = stephist(incidents[incidents .< medi + spread*2], bins=0:1: medi + spread*2, color=:black, label=:none, alpha=1, size=(500,300), dpi=1000)
#histogram!(fiddiincidents[fiddiincidents .<  medi + spread*3], bins = 0:1: medi + spread*3, color=:blue, label="50ps", alpha=0.3)
#histogram!(hunniincidents[hunniincidents .<  medi + spread*3], bins = 0:1: medi + spread*3, color=:red, label="100ps", alpha=0.3)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")

cutread = CSV.read("../cutparams.csv", DataFrame; header=false, delim=",", ignorerepeated=false)
cutparams = collect(cutread[1,:])

#cut(E, p) = p[1] .+ p[2] * exp.(- (E .+ p[3]) ./ p[4])
cut(E, p) = p[1] .+ p[2] ./ (E .+ p[3])
cutlee = 0.5 # The cut isn't perfect, so add a bit of extra leeway
lowert = 0.5 # Set a lower threshold for detection. Usefull for cutting out crosstalk also

cutindices = cutindices = firsts .< cut(incidents, cutparams) .+ cutlee .&& firsts .> lowert .&& seconds .> 0.
cutincidents = incidents[cutindices]

extracutindices = firsts .< cut(incidents, cutparams) .+ cutlee .&& firsts .> lowert .&& seconds .> lowert .&&  seconds .< cut(incidents, cutparams) .+ cutlee
extracutincidents = incidents[extracutindices]

histogram!(cutincidents[cutincidents .<  medi + spread*2],  bins=0:1: medi + spread*2, alpha=0.9, color=:red, label="Cut")
histogram!(extracutincidents[extracutincidents .<  medi + spread*2],  bins=0:1: medi + spread*2, alpha=0.9, color=:blue, label="Extra cut")
savefig("plots/TotalEnergiesCut.svg")
display(fig2)




fig3 = histogram2d(incidents[incidents .< medi*2], firsts[incidents .< medi*2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")

fitpoints = collect(25:maximum(incidents[incidents .< medi *2]))
plot!(fitpoints, cut(fitpoints, cutparams) .+ cutlee, label="Upper cut")
hline!([lowert], label="Lower cut")
savefig("plots/FirstHeatmapCut.svg")
#display(fig3)

fig4 = histogram2d(incidents[incidents .< medi*2], seconds[incidents .< medi*2], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.svg")
#display(fig4)

# Writing into results file
resultfile = open("results.csv", "a")
if filesize("results.csv") == 0 # Checks if the file is empty, and writes a header if it is
    write(resultfile, "Median, spread, 50ps median, 50ps spread, 100ps median, 100ps spread\n")
end
write(resultfile, string(medi),", ", string(spread),", ", string(fiddimedi),", ", string(fiddispread),", ", string(hunnimedi),", ", string(hunnispread), "\n")
close(resultfile)



function signalnoise(inci, signalstart, signalstop)
	indices = inci .> signalstart .&& inci .< signalstop
	
	return [sum(indices), sum(.! indices), sum(indices)/ sum(.! indices) ]
end


println("No cut: #",length(incidents),", S,N,S/N: ", signalnoise(incidents, 95,105))

println("Cut: #",length(cutincidents),", S,N,S/N: ", signalnoise(cutincidents, 95,105))

println("Extra cut: #",length(extracutincidents),", S,N,S/N: ", signalnoise(extracutincidents, 95,105))
