
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

using LsqFit

# Reading data
#data = data = CSV.read("100milSWMinus_PhaseSpace.phsp", DataFrame; header=false, delim=" ", ignorerepeated=true)
data = data = CSV.read("SWMinus_PhaseSpace.phsp", DataFrame; header=false, delim=" ", ignorerepeated=true)

#println(data)


# Filtering for neutrons and finding the energy distribution
neutrons = data[:,8] .== 2112
bins = range(0, maximum(data[neutrons, 6]), length=250)
fig0 = histogram(data[neutrons,6], bins=bins)
title!("Spallation neutron energies")
xlabel!("Energy (MeV)")
ylabel!("Counts in bin")
savefig("NeutronEnergies.svg")
display(fig0)

#println(neutrons)
# Filtering for gammas and finding the energy distribution
gammas = data[:,8] .== 22
bins = range(0, maximum(data[gammas, 6]), length=250)
fig1 = histogram(data[gammas,6], bins=bins)
title!("Gamma energies")
xlabel!("Energy (MeV)")
ylabel!("Counts in bin")
savefig("GammaEnergies.svg")
display(fig1)

others = data[:,8] .== 22 .&& data[:,8] .== 2112


# Positions of the neutrons on 15cm x 15cm square at 1m distance roughly
fig2 = histogram2d(data[neutrons,1], data[neutrons,2], bins=(50,50))
title!("Neutron XY heatmap")
xlabel!("X (cm)")
ylabel!("Y (cm)")
savefig("NeutronHeatmap.svg")
display(fig2)


# X and Y histogram, to check whether they are distributed uniformly or not
fig3 = stephist(data[neutrons,1], bins=100,label="X")
stephist!(data[neutrons,2], bins=100,label="Y")
title!("X and Y histograms")
xlabel!("Position (cm)")
ylabel!("Counts in bin")
savefig("NeutronXYHistograms.svg")
display(fig3)


# Checking the spatial distribution vs the direction (they should preferably all point back to roughly the same point)
energindex = data[:,6] .> 1
ds = sqrt.(data[neutrons,1].^2 .+ data[neutrons,2].^2)
ts = sqrt.(data[neutrons,4].^2 .+ data[neutrons,5].^2)

fig4 = histogram2d(ts[ts .< 0.3], ds[ts .< 0.3], bins=(500,500), label="Simulation")
title!("Angular distribution of neutrons")
ylabel!("Distance from center (cm)")
xlabel!("Momentum direction away from center (radians)")
xlims!(0, 0.3)
ylims!(0, 22)
model(x, p) = sin.(x) .* p[1]
fit = curve_fit(model, ts, ds, [100.])
plot!(x->sin(x) .* 100, label="sin(θ)⋅100cm", color=:red, linestyle=:dash)
#plot!(x->sin(x) .* coef(fit)[1], label="fit to sin(θ)⋅A", color=:red)
savefig("NeutronDirection.svg")
display(fig4)


fig5 = histogram2d(data[neutrons, 6][ts .< 0.3], ts[ts .< 0.3], bins=(200,200))
title!("Direction and energy of neutrons")
xlabel!("Energy (MeV)")
ylabel!("Momentum direction away from center (radians)")
savefig("DirectionEnergy.svg")
display(fig5)

fig6 = histogram2d(data[neutrons,6][ts .< 0.3], ts[ts .< 0.3] ./ ds[ts .< 0.3], bins=(200,200))
title!("Non-isotropicness")
savefig("Enertropic.svg")





println(length(data[neutrons,1]), " neutrons, ", length(data[gammas,1]), " gammas.")
nx = data[neutrons,1]
ny = data[neutrons,2]

centisquare = abs.(data[neutrons,1]) .< 1 .&& abs.(data[neutrons,2]) .< 1
println(sum(centisquare))
