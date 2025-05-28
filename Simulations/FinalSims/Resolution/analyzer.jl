
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
medi2=data[:,2]
spread=data[:,3]
spreads=data[:,4]
thirty=data[:,5]
thirtys=data[:,6]
thirtyspread=data[:,7]
thirtyspreads=data[:,8]
fiddi=data[:,9]
fiddis=data[:,10]
fiddispread=data[:,11]
fiddispreads=data[:,12]
hunni=data[:,13]
hunnis=data[:,14]
hunnispread=data[:,15]
hunnispreads=data[:,16]

println(medi)


fig1 = plot(energies,100* spread ./ energies, yerr = 100* spreads ./ energies, c=:black, markerstrokecolor=:black, label="Ideal")
plot!(energies,100 * thirtyspread ./ energies, yerr = 100* thirtyspreads ./ energies, c=:green, markerstrokecolor=:green, label="30ps")
plot!(energies,100 * fiddispread ./ energies, yerr = 100* fiddispreads ./ energies, c=:blue, markerstrokecolor=:blue, label="50ps")
plot!(energies,100 * hunnispread ./ energies, yerr = 100* hunnispreads ./ energies, c=:red, markerstrokecolor=:red, label="100ps")


title!("Uncertainty of energy calculation")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy uncertainty (%)")
savefig("plots/uncertainty.png")
#display(fig1)


