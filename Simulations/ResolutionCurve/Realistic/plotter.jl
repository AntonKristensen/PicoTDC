
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

fiddidata = CSV.read("output/fiddiuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)
fiddiincidents = fiddidata[:,1]

hunnidata = CSV.read("output/hunniuncertainmatches.csv", DataFrame; header=1, delim=",", ignorerepeated=false)
hunniincidents = hunnidata[:,1]


###############
# Collecting results from all individual pairs into one big list
incidents = data[:,1]
firsts = data[:,2]
seconds = data[:,3]
fronts = data[:,4]
backs = data[:,5]


function statisticing(data)
    # Doing a slight bit of statistics
    medi = median(data) # Getting a decent robust guess for the mean of the peak so I can make a un-bad cut when fitting
    lower = medi - quantile(data, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
    upper = quantile(data, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak

    #print(medi, ", ", lower, ", ", upper, "\n")

    bound = min(lower,upper)

    # ML fit, cutting data 2 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound*2) .& (data .< medi + bound*2)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    print(medi, " ", bound, " ", params(gaussfit), "\n")

    # ML fit, cutting data 1 sigma below and above the calculated mean
    fitdata = data[(data .> medi - bound*2) .& (data .< medi + bound*2)] # Cutting a roughly 3sigma region around the peak
    gaussfit = fit_mle(Normal, fitdata)
    print(medi, ", ", lower, ", ", upper, " ", params(gaussfit),"\n")

    return params(gaussfit)
end

medi, spread = statisticing(incidents)
fiddimedi, fiddispread = statisticing(fiddiincidents)
hunnimedi, hunnispread = statisticing(hunniincidents)



fig2 = histogram(incidents[incidents .< medi*2], bins = 250, color=:black, label="Ideal", alpha=1, size=(400,600))
histogram!(fiddiincidents[fiddiincidents .< medi*2], bins = 250, color=:blue, label="50ps", alpha=0.3)
histogram!(hunniincidents[hunniincidents .< medi*2], bins = 250, color=:red, label="100ps", alpha=0.3)


title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.png")
savefig("plots/TotalEnergies.pdf")
#display(fig2)


fig3 = histogram2d(incidents[incidents .< medi*2], firsts[incidents .< medi*2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.png")
#display(fig3)

fig4 = histogram2d(incidents[incidents .< medi*2], seconds[incidents .< medi*2], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)

# Writing into results file
resultfile = open("results.csv", "a")
if filesize("results.csv") == 0 # Checks if the file is empty, and writes a header if it is
    write(resultfile, "Median, spread, 50ps median, 50ps spread, 100ps median, 100ps spread\n")
end
write(resultfile, string(medi),", ", string(spread),", ", string(fiddimedi),", ", string(fiddispread),", ", string(hunnimedi),", ", string(hunnispread), "\n")
close(resultfile)
