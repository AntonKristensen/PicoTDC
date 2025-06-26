
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using LsqFit


filepath = "stochasticresults.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

sort!(data, names(data)[1])
println(data)

flux = 1 ./ (1.0e-12 * data[:,1] * 6^2)

fig1 = plot(flux, data[:,4], linewidth=3, legend=:none,yscale = :log10,xscale = :log10, minorgrid = true)
title!("Signal/Noise above 25 MeV")
xlabel!("Neutron flux (s⁻¹ cm⁻²)")
ylabel!("S/N")
display(fig1)
savefig("stochasticnoise.svg")

