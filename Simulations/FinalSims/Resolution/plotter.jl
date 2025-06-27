
using Plots
using CSV
using DataFrames
using Statistics
using StatsBase
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

thirtydata = CSV.read("output/thirtyuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

fiddidata = CSV.read("output/fiddiuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

hunnidata = CSV.read("output/hunniuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)

include("../statistic.jl")

include("../cutting.jl")
cutdata = cut(data)
cutthirtydata = cut(thirtydata)
cutfiddidata = cut(fiddidata)
cuthunnidata = cut(hunnidata)


###############
# Collecting results from all individual pairs into one big list

med, spread, amount, region = bootstatisticing(cutdata[:,1], 10000)
thirtymedi, thirtyspread, thirtyamount, thirtyregion = bootstatisticing(cutthirtydata[:,1], 10000)
fiddimedi, fiddispread, fiddiamount, fiddiregion = bootstatisticing(cutfiddidata[:,1], 10000)
hunnimedi, hunnispread, hunniamount, hunniregion = bootstatisticing(cuthunnidata[:,1], 10000)

range = 1.2
medi = med[1]
i = cutdata[:,1] .< medi*range
points = collect(region[1] - region[2]: 0.1 : region[1] + region[2])
thirtypoints = collect(thirtyregion[1] - thirtyregion[2]: 0.1 : thirtyregion[1] + thirtyregion[2])
fiddipoints = collect(fiddiregion[1] - fiddiregion[2]: 0.1 : fiddiregion[1] + fiddiregion[2])
hunnipoints = collect(hunniregion[1] - hunniregion[2]: 0.1 : hunniregion[1] + hunniregion[2])


println(amount)

fig2 = stephist(cutdata[:,1], bins=0:1:medi*range, color=:black, label="Ideal", alpha=1, size=(500,300), dpi=1000)
plot!(points, amount[1] *  pdf(Normal(med[1], spread[1]), points), label="Ideal fit", color=:black)
histogram!(cutthirtydata[:,1], bins=0:1:medi*range, color=:green, label="30ps", alpha=0.50)
plot!(thirtypoints, thirtyamount[1] *  pdf(Normal(thirtymedi[1], thirtyspread[1]), thirtypoints), label="30ps fit", color=:green)
histogram!(cutfiddidata[:,1], bins=0:1:medi*range, color=:blue, label="50ps", alpha=0.50)
plot!(fiddipoints, fiddiamount[1] *  pdf(Normal(fiddimedi[1], fiddispread[1]), fiddipoints), label="50ps fit", color=:blue)
histogram!(cuthunnidata[:,1], bins=0:1:medi*range, color=:red, label="100ps", alpha=0.50)
plot!(hunnipoints, hunniamount[1] *  pdf(Normal(hunnimedi[1], hunnispread[1]), hunnipoints), label="100ps fit", color=:red)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
display(fig2)

fig2zoom = stephist(cutdata[:,1], bins=hunnipoints, color=:black, label="Ideal", alpha=1, size=(500,300), dpi=1000)
plot!(points, 0.1 *  amount[1] *  pdf(Normal(med[1], spread[1]), points), label="Ideal fit", color=:black, linewidth=3)
histogram!(cutthirtydata[:,1], bins=hunnipoints, color=:green, label="30ps", alpha=0.50)
plot!(thirtypoints,  0.1 * thirtyamount[1] *  pdf(Normal(thirtymedi[1], thirtyspread[1]), thirtypoints), label="30ps fit", color=:green, linewidth=3)
histogram!(cutfiddidata[:,1], bins=hunnipoints, color=:blue, label="50ps", alpha=0.50)
plot!(fiddipoints,  0.1 * fiddiamount[1] *  pdf(Normal(fiddimedi[1], fiddispread[1]), fiddipoints), label="50ps fit", color=:blue, linewidth=3)
histogram!(cuthunnidata[:,1], bins=hunnipoints, color=:red, label="100ps", alpha=0.50)
plot!(hunnipoints,  0.1 * hunniamount[1] *  pdf(Normal(hunnimedi[1], hunnispread[1]), hunnipoints), label="100ps fit", color=:red, linewidth=3)
title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergiesZoomed.svg")
display(fig2zoom)




fig3 = histogram2d(cutdata[i,1], cutdata[i,2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.png")
#display(fig3)

fig4 = histogram2d(cutdata[i,1], cutdata[i,3], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)

# Writing into results file
resultfile = open("results.csv", "a")
if filesize("results.csv") == 0 # Checks if the file is empty, and writes a header if it is
    write(resultfile, "Median, sigma, spread, sigma, 30ps median, sigma, 30ps spread, sigma, 50ps median, sigma, 50ps spread, sigma, 100ps median, sigma, 100ps spread, sigma\n")
end
write(resultfile, string(med[1]),", ", string(med[2]),", ", string(spread[1]),", ", string(spread[2]),", ", string(thirtymedi[1]),", ", string(thirtymedi[2]),", ", string(thirtyspread[1]),", ", string(thirtyspread[2]),", ", string(fiddimedi[1]),", ", string(fiddimedi[2]),", ", string(fiddispread[1]),", ", string(fiddispread[2]),", ", string(hunnimedi[1]),", ", string(hunnimedi[2]),", ", string(hunnispread[1]),", ", string(hunnispread[2]), "\n")
close(resultfile)

