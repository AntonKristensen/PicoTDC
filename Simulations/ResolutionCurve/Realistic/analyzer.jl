
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
fiddi=data[:,3]
fiddispread=data[:,4]
hunni=data[:,5]
hunnispread=data[:,6]
two=data[:,7]
twospread=data[:,8]
tree=data[:,9]
treespread=data[:,10]



fig1 = plot(energies, spread ./ energies, c=:black, label="Ideal")
#plot!(energies, energies - medi, c=:black, label="Syst Ideal")
plot!(energies,fiddispread ./ energies, c=:blue, label="50ps")
#plot!(energies, energies - fiddi, c=:blue, label="syst 50ps")
plot!(energies,hunnispread ./ energies, c=:red, label="100ps")
#plot!(energies, energies - hunni, c=:red, label="Syst 100ps")
plot!(energies,twospread ./ energies, c=:green, label="200ps")
#plot!(energies, energies - two, c=:green, label="Syst 200ps")
plot!(energies,treespread ./ energies, c=:magenta, label="350ps")
#plot!(energies, energies - tree, c=:magenta, label="Syst 350ps")
title!("Energy uncertainty")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Uncertainty in energy measurement (MeV)")
savefig("plots/uncertainty.png")
#display(fig1)


