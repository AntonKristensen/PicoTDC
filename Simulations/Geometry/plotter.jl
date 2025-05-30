
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)



###############
# Collecting results from all individual pairs into one big list
incidents = data[:,1]
firsts = data[:,2]
seconds = data[:,3]
fronts = data[:,4]
backs = data[:,5]
frontevents = data[:,6]
backevents = data[:,7]

cincidents = incidents[frontevents .== backevents]
fincidents = incidents[frontevents .!= backevents]

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
    return params(gaussfit)
end

medi, spread = statisticing(incidents)


fig2 = histogram(cincidents, bins = 0:1:300, color=:blue, label="Correct events", alpha=0.5)
histogram!(fincidents, bins = 0:1:300, color=:red, label="False events", alpha=0.5)
title!("Total energy spectrum, 1.5e5 n/cm²/s")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.png")
#display(fig2)


fig3 = histogram2d(incidents[incidents .< medi*2], firsts[incidents .< medi*2], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")
#display(fig3)

fig4 = histogram2d(incidents[incidents .< 300], seconds[incidents .< 300], bins=(250, 250))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.png")
#display(fig4)

println("Correct: ", length(incidents[frontevents .== backevents]), ", False: ", length(incidents[frontevents .!= backevents]))
