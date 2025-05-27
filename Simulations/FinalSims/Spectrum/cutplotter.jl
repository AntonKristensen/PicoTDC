
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using LsqFit


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





fig4 = histogram2d(incidents[incidents .< 250], seconds[incidents .< 250], bins=(150, 150))
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.svg")

x = [28, 38, 50, 75, 100, 125, 150, 175, 200, 225]
y = [23, 14, 9, 6, 4.7, 4, 3.2, 3, 2.85, 2.75]

alt(E, p) = p[1] .+ p[2] ./ (E .+ p[3])
altfit = curve_fit(alt, x, y, [5., 100., -5.])

cut(E, p) = p[1] .+ p[2] * exp.(- (E .+ p[3]) ./ p[4])
p0 = [5.0, 25.0, 1.0, 100.0]
expfit = curve_fit(cut, x, y, p0)
fitpoints = collect(25:250)
scatter!(x,y, color=:green, alpha=0.8, label="")
#plot!(fitpoints, cut(fitpoints, expfit.param), label="Exponential fit", color=:green) 
plot!(fitpoints, alt(fitpoints, altfit.param), label="Hyperboloid fit", color=:green, alpha=0.8, linewidth=2)
savefig("plots/SecondHeatmapFit.svg")
display(fig4)




fig3 = histogram2d(incidents[incidents .< 250], firsts[incidents .< 250], bins=(150, 150))
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")

plot!(fitpoints, cut(fitpoints, expfit.param) .+ 1, label="Upper cut", color=:orange, alpha=0.8, linewidth=2)
hline!([0.5], label="Lower cut", color=:orange, alpha=0.8, linewidth=2)
savefig("plots/FirstHeatmapCut.svg")
#display(fig3)



fig2 = histogram(incidents[incidents .< 250], bins = 250, color=:black, label="No cut", alpha=1, size=(500,300), dpi=1000)
title!("Neutron Spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("plots/TotalEnergies.svg")

uppercutindices = firsts .< cut(incidents, expfit.param) .+ 1
uppercutincidents = incidents[uppercutindices]
histogram!(uppercutincidents[uppercutincidents .< 250], bins=250, label="Upper cut", color=:red, alpha=0.9)

cutindices = uppercutindices .&& firsts .> 0.5 .&& seconds .> 0.5
cutincidents = incidents[cutindices]
histogram!(cutincidents[cutincidents .< 250], bins=250, label="Upper+lower cut", color=:blue, alpha=0.9)

savefig("plots/TotalEnergiesCut.svg")
#savefig("plots/TotalEnergies.pdf")
#display(fig2)



function savecut(cutparams)
	file = open("../cutparams.csv", "w")
	for n in 1:length(cutparams)
		if n == length(cutparams)
			write(file, string(cutparams[n]))
		else
			write(file, string(cutparams[n]) * ", ")
		end
	end
	close(file)
end

savecut(altfit.param)





