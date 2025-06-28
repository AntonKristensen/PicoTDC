
using Plots
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions
using Glob

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

include("../statistic.jl")

include("../cutting.jl")
#cutdata = cut(data)
cutdata = data

###############
# Collecting results from all individual pairs into one big list



#medi, spread = statisticing(cutdata[:,1])



phasefiles = glob("output/PhaseSpaces/*.phsp")

energies = Float64[]
particles = Int[]
for file in phasefiles
    println(file)
    phasespaces = CSV.read(file, DataFrame; delim=" ", ignorerepeated=true, ignoreemptyrows=true)
    e = phasespaces[:,6]
    append!(energies, e)
    t = phasespaces[:,8]
    append!(particles, t)
end

neutrons = particles .== 2112
protons = particles .== 2212
electrons = particles .== 11
gammas = particles .== 22


fig1 = stephist(energies[neutrons], bins=1:1:250, color=:black, label="Neutrons: "*string(sum(neutrons)))
stephist!(energies[protons], bins=1:1:250, color=:red, label="Protonss: "*string(sum(protons)))
stephist!(energies[electrons], bins=1:1:250, color=:blue, label="Electronss: "*string(sum(electrons)))
stephist!(energies[gammas], bins=1:1:250, color=:green, label="Gammas: "*string(sum(gammas)))
xlabel!("Energy (MeV)")
ylabel!("Counts")
title!("Particles at detector")
savefig("plots/Energies.svg")
display(fig1)



fig2 = stephist(cutdata[:,1], bins=1:1:250, color=:black, label="Data", alpha=1, size=(500,300), dpi=1000)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
display(fig2)




fig3 = histogram2d(cutdata[:,1], cutdata[:,2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.png")
#display(fig3)

fig4 = histogram2d(cutdata[:,1], cutdata[:,3], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)



