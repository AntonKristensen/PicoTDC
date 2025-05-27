
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



fig2 = histogram(incidents[incidents .< 250], bins = 250, color=:black, label="Ideal", alpha=1, size=(500,300), dpi=1000)


title!("Monoenergetic Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
#savefig("plots/TotalEnergies.pdf")
#display(fig2)


fig3 = histogram2d(incidents[incidents .< 250], firsts[incidents .< 250], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")
#display(fig3)

fig4 = histogram2d(incidents[incidents .< 250], seconds[incidents .< 250], bins=(150, 150))

x = [28, 38, 50, 75, 100, 125, 150]
y = [20, 15, 9, 6, 5.5, 4, 3.2]

scatter!(x,y, color=:green, alpha=0.7)
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.svg")
#display(fig4)



