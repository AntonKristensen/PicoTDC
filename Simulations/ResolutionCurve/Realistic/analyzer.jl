
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "results.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

start = parse(Int,ARGS[1])
interval = parse(Int,ARGS[2])
stop = parse(Int,ARGS[3])

println(start, ", ", interval, ", ", stop)

energies = collect(start:interval:stop)

println(energies)

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

title!("Uncertainty of energy calculation")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy uncertainty (%)")
savefig("plots/uncertainty.png")
#display(fig1)


