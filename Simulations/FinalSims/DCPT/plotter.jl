
using Plots
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

include("../statistic.jl")

include("../cutting.jl")
cutdata = cut(data)


###############
# Collecting results from all individual pairs into one big list



#medi, spread = statisticing(cutdata[:,1])



fig2 = stephist(cutdata[:,1], bins=1:1:250, color=:black, label="Data", alpha=1, size=(500,300), dpi=1000)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
display(fig2)




fig3 = histogram2d(cutdata, cutdata, bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.png")
#display(fig3)

fig4 = histogram2d(cutdata, cutdata, bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)

