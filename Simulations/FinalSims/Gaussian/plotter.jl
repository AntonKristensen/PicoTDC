
using Plots
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

#thirtydata = CSV.read("output/thirtyuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

fiddidata = CSV.read("output/fiddiuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

#hunnidata = CSV.read("output/hunniuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

include("../statistic.jl")

include("../cutting.jl")
cutdata = cut(data)
#cutthirtydata = cut(thirtydata)
cutfiddidata = cut(fiddidata)
#cuthunnidata = cut(hunnidata)


###############
# Collecting results from all individual pairs into one big list



media, spread = bootstatisticing(cutdata[:,1])
#thirtymedi, thirtyspread = statisticing(cutthirtydata[:,1])
fiddimedi, fiddispread = statisticing(cutfiddidata[:,1])
#hunnimedi, hunnispread = statisticing(cuthunnidata[:,1])
medi = media[1]
range = 1.5
i = cutdata[:,1] .< medi*range
normalarea = length(cutdata[cutdata[:,1] .< 140 .&& cutdata[:,1] .> 60,1])

fig2 = histogram(cutdata[:,1], bins=0:1:medi*range, color=:black, label="Ideal", alpha=1, size=(500,300), dpi=1000, legend=:topleft)
#histogram!(cutthirtydata[:,1], bins=0:1:medi*range, color=:green, label="30ps", alpha=0.5)
#histogram!(cutfiddidata[:,1], bins=0:1:medi*range, color=:blue, label="50ps", alpha=0.5)
#histogram!(cuthunnidata[:,1], bins=0:1:medi*range, color=:red, label="100ps", alpha=0.5)
plot!(60:1:140, normalarea * pdf(Normal(100,20), 60:1:140), label="Normal, μ=100MeV, σ=20MeV)", linecolor=:red, linewidth=3)
title!("Gaussian Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
display(fig2)




fig3 = histogram2d(cutdata[i,1], cutdata[i,2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.png")
#display(fig3)

fig4 = histogram2d(cutdata[i,1], cutdata[i,3], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)

#=
# Writing into results file
resultfile = open("results.csv", "a")
if filesize("results.csv") == 0 # Checks if the file is empty, and writes a header if it is
    write(resultfile, "Median, spread, 30ps median, 30ps spread, 50ps median, 50ps spread, 100ps median, 100ps spread\n")
end
write(resultfile, string(medi),", ", string(spread),", ", string(thirtymedi),", ", string(thirtyspread),", ", string(fiddimedi),", ", string(fiddispread),", ", string(hunnimedi),", ", string(hunnispread), "\n")
close(resultfile)
=#
