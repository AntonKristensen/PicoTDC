
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
thirty=data[:,3]
thirtyspread=data[:,4]
fiddi=data[:,5]
fiddispread=data[:,6]
hunni=data[:,7]
hunnispread=data[:,8]



fig1 = plot(energies,100* spread ./ energies, c=:black, label="Ideal")
plot!(energies,100 * thirtyspread ./ energies, c=:green, label="30ps")
plot!(energies,100 * fiddispread ./ energies, c=:blue, label="50ps")
plot!(energies,100 * hunnispread ./ energies, c=:red, label="100ps")


title!("Uncertainty of energy calculation")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy uncertainty (%)")
savefig("plots/uncertainty.png")
#display(fig1)


