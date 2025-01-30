
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "results.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

energies = collect(10:10:250)

##############
medi=data[:,1]
spread=data[:,2]

fig1 = plot(energies, spread)
title!("Energy uncertainty")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Uncertainty in energy measurement (MeV)")
savefig("plots/uncertainty.png")
#display(fig1)

