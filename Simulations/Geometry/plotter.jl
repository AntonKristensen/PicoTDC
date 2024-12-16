
using Plots
using CSV
using DataFrames
using Statistics
using Distributions


data = CSV.read("output/results.csv", DataFrame; delim=",", ignorerepeated=true)

energies = unique(data[:,1]) # Extract the various energies simulated
means = zeros(length(energies))
meanerrors = zeros(length(energies))
stds = zeros(length(energies))
stderrors = zeros(length(energies))

for i in eachindex(energies)
   means[i] = mean(data[data[:,1] .== energies[i],2])
   meanerrors[i] = std(data[data[:,1] .== energies[i],2]) / sqrt(length(data[data[:,1] .== energies[i],2]))
   stds[i] = mean(data[data[:,1] .== energies[i],3])
   stderrors[i] = std(data[data[:,1] .== energies[i],3]) / sqrt(length(data[data[:,1] .== energies[i],3]))
end

fig1 = plot(energies, energies .- means, yerr = meanerrors, markershape=:circle, linetype=:scatter, label="Simulation")
plot!(energies, energies .- means, ribbon = meanerrors, fillalpha = 0.15, c = 1, linewidth=0, label = "Confidence band")

title!("Energy bias")
xlabel!("Energy of incident neutrons (MeV)")
ylabel!("Bias in determined energy (MeV)")
savefig("plots/EnergyBias.svg")
display(fig1)

fig2 = plot(energies, stds, yerr = stderrors, markershape=:circle, linetype=:scatter, label="Simulation")
plot!(energies, stds, ribbon = stderrors, fillalpha = 0.15, c = 1, linewidth=0, label = "Confidence band")
title!("Energy uncertainty")
xlabel!("Energy of incident neutrons (MeV)")
ylabel!("Uncertainty in determined energy (MeV)")
savefig("plots/EnergyUncertainty.svg")
display(fig2)
