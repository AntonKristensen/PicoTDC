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



fig1 = plot(scintillators * 2, correct ./ flux, linewidth=3, legend=:false) 
title!("Efficiency")
ylabel!("Matches per neutron flux (cm²)")
xlabel!("Distance between arrays (cm)")
savefig("plots/distance.svg")
display(fig1)

