
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "results.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

energies = collect(25:25:250)

##############
medi=data[:,1]
spread=data[:,2]
fiddi=data[:,3]
fiddispread=data[:,4]
hunni=data[:,5]
hunnispread=data[:,6]



fig1 = plot(energies,100* spread ./ energies, c=:black, label="Ideal")
#plot!(energies, energies - medi, c=:black, label="Syst Ideal")
plot!(energies,100 * fiddispread ./ energies, c=:blue, label="50ps")
#plot!(energies, energies - fiddi, c=:blue, label="syst 50ps")
plot!(energies,100 * hunnispread ./ energies, c=:red, label="100ps")
#plot!(energies, energies - hunni, c=:red, label="Syst 100ps")

title!("Energy uncertainty")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy uncertainty (%)")
savefig("plots/uncertainty.png")
#display(fig1)


