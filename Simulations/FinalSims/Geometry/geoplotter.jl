using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using LsqFit


filepath = "geometryresults.txt"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

scintillators = data[:,1]
correct = data[:,2]
fake = data[:,3]
neutrons = data[:,4]
flux = data[:,5]

println(neutrons ./ flux)


fig1 = plot(scintillators, correct  ./ flux ./ 1 ) 
title!("Efficiency")
ylabel!("Matches per neutron flux (match cm²)")
xlabel!("Scintillator pairs")
savefig("plots/geometry.svg")
display(fig1)

fig2 = plot(scintillators, correct  ./ flux ./ (scintillators * 0.6^2) ) 
title!("Efficiency per scintillator pair")
ylabel!("Area normalized matches")
xlabel!("Scintillator pairs")
savefig("plots/normalizedgeometry.svg")
display(fig2)