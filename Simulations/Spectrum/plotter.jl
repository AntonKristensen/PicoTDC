
using Plots
using CSV
using DataFrames
using Statistics
using Distributions

filepath = "output/matches.csv"
data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

print(data)

###############
# Collecting results from all individual pairs into one big list
incidents = data[:,1]
firsts = data[:,2]
seconds = data[:,3]
fronts = data[:,4]
backs = data[:,5]


# Doing a slight bit of statistics
medi = median(incidents) # Getting a decent robust guess for the mean of the peak so I can make a un-bad cut when fitting
lower = medi - quantile(incidents, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
upper = quantile(incidents, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak

# ML fit, cutting data 20% below and above the calculated mean
fitdata = incidents[(incidents .> medi - lower*3) .& (incidents .< medi + upper*3)] # Cutting a roughly 3sigma region around the peak
gaussfit = fit_mle(Normal, fitdata)
print(params(gaussfit), (upper + lower)/2, "\n")

#fig1 = plot()
#for i in 1:length(incidentframe[:,3])
#    if length(incidentframe[i,3]) > 0 # It breaks if it there weren't any matches
#        stephist!(incidentframe[i,3], label="Front"*string(incidentframe[i,1])*", Back"*string(incidentframe[i,2]))
#        #print(length(incidentframe[i,3]),"\n")
#    end
#end
#
#xlims!(0,medi*2)
#savefig("plots/SeparateEnergies.svg")
#title!("Energies from detector pairs")
#xlabel!("Energy (MeV)")
#ylabel!("Counts")
#savefig("plots/SeparateEnergies.svg")
#display(fig1)

fig2 = histogram(incidents[incidents .< medi*2], bins = 100)
title!("Total energy spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")
display(fig2)


fig3 = histogram2d(incidents[incidents .< medi*2], firsts[incidents .< medi*2], bins=(50, 50))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")
display(fig3)

fig4 = histogram2d(incidents[incidents .< medi*2], seconds[incidents .< medi*2], bins=(50, 50))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.svg")
display(fig4)

# Writing into results file
file = open("output/results.csv", "a")
write(file, string(params(gaussfit)[1]) * ", " * string(params(gaussfit)[2]))
close(file)
