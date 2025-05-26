

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


include("functions.jl")

fil = "output/front1.phsp"

readed = reader(fil)

#@time collected = collector(fil)
#print(collected)

elast = readed[:,4] .== "hadElastic"
inelast = readed[:,4] .== "neutronInelastic"
prot = readed[:,2] .== "proton"
neut = readed[:,2] .== "neutron"
deut = readed[:,2] .== "deuteron" .|| readed[:,2] .== "triton"
exits = readed[:,8] .> 0

fig1 = stephist(readed[elast .&& prot .&& exits, 8], bins=100, label="Elastic protons", alpha=0.75, linecolor=:black, linewidth=2)
stephist!(readed[inelast .&& exits, 8], bins=100, label="All Fragments", alpha=0.75, linecolor=:green, linewidth=2)
stephist!(readed[inelast .&& exits .&& prot, 8], bins=100, label="Proton fragments", alpha=0.75, linecolor=:red, linewidth=2)
stephist!(readed[inelast .&& exits .&& neut, 8], bins=100, label="Neutron fragments", alpha=0.75, linecolor=:blue, linewidth=2)
stephist!(readed[inelast .&& exits .&& deut, 8], bins=100, label="Deuteron & triton fragments", alpha=0.75, linecolor=:magenta, linewidth=2)
stephist!(readed[inelast .&& exits .&& .!neut .&& .!prot .&& .!deut, 8], bins=100, label="Nuclear fragments", alpha=0.75, linecolor=:orange, linewidth=2)
title!("Particle exit energies 200 MeV")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/InelasticExitEnergies.svg")
display(fig1)

fig2 = stephist(readed[elast .&& prot .&& exits, 8], bins=100, label="Elastic protons", alpha=0.75, linecolor=:black, linewidth=2)
stephist!(readed[elast .&& exits .&& neut, 8], bins=100, label="Neutron fragments", alpha=0.75, linecolor=:blue, linewidth=2)
title!("Elastic exit energies")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/ElasticExitEnergies.svg")
#display(fig2)



#print(readed[readed[:,8] .> 90 .&& inelast,:])
#println(readed[readed[:,8] .> 50 .&& inelast .&& prot])

#println(readed[readed[:,3] .== 49977014, :])
#print(reader(fil))
#print(collected)
#incidentframe = detectorlooping(minimumenergy=1, addedtime = 1e-6)
#matchwriter(incidentframe)







