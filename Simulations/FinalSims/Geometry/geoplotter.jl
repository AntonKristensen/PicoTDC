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


fig1 = plot(scintillators, correct .* neutrons ./ (flux .* scintillators))
title!("Efficiency per scintillator")
ylabel!("Matches per neutron per scintillator pair")
xlabel!("Scintillator pairs")
savefig("plots/geometry.svg")
display(fig1)

