
using Plots
using CSV
using DataFrames
using Statistics
using Distributions


# Making a function to collect energy deposit from same event but different particles (such as electron secondaries or inelastic products)
function collector(filepath)
    # Reading file with data in it
    data = CSV.read(filepath, DataFrame; header=false, delim=" ", ignorerepeated=true)

    # Making a list of indeces to be popped, and adding the extra energy deposit from the extra particles (usually electrons)
    # Also needs to make sure that the time stamp is the lowest of all the created particles, because Geant4 does not for some reason spit them out in order
    pops = []
    for i in length(data[:,1]):-1:1+1
        if data[i, 3] == data[i-1, 3] # If the particles are in the same Geant4Run
            append!(pops, i)
            data[i-1, 1] += data[i,1] # Add the energy together
            data[i-1, 9] = minimum([data[i-1, 9], data[i, 9]]) # Choose the lowest time stamp
        end
    end

    delete!(data, reverse(pops)) # Deleting the secondaries
    return data
end


# Plotting some histograms
function histogrammer(filepath)
    data = collector(filepath)
    bins = range(0, maximum(data[:, 1]), length=100)

    fig = histogram(data[:, 1], bins = bins, label="Total", normalize=:false, color=:black, alpha=0.2)
    title!("Energy deposition in first scintillator")
    xlabel!("E (MeV)")
    ylabel!("Counts")

    protons = filter(row -> row.Column2 == "proton", filter(row -> row.Column4 == "hadElastic", data))
    stephist!(protons[:, 1], bins = bins, label="Elastic Protons", color=:red)

    inelastics = filter(row -> row.Column4 == "neutronInelastic", data)
    stephist!(inelastics[:, 1], bins = bins, label="Inelastics", color=:green)

    carbons = filter(row -> row.Column2 == "C12", data)
    stephist!(carbons[:, 1], bins = bins, label="C12", color=:blue)

    display(fig)

    plotfilename = filepath[1 : findlast(".", filepath)[1]] * "svg"
    savefig(plotfilename)
end

#histogrammer("First.phsp")
#histogrammer("Second.phsp")

# Adding time, some nanoseconds between each beam neutron
addedtime = 1e-9

first = collector("First.phsp")
first[!,"Column10"] = first[:,3]*addedtime + first[:,end]

second = collector("Second.phsp")
second[!,"Column10"] = second[:,3]*addedtime + second[:,end]

sort!(first, [:Column10])
sort!(second, [:Column10])


# Function that matches hits within some time difference
function matcher(a, b, mindif, maxdif)
    aindex = 1
    matches = []
    for bindex in 1:length(b[:, end])

        # Finding the first index where the time difference is more than the minimum
        while b[bindex, end] - a[aindex, end] > maxdif
            aindex += 1
        end        

        upperindex = aindex 
        while b[bindex, end] - a[upperindex, end] > mindif
            push!(matches, [upperindex, bindex])
            upperindex += 1
        end

        #print(aindex, "\n")
        #print(b[bindex, end])
    end
    return matches
end

matches = matcher(first, second, 1.7e-9, 10e-9) # 10ns max time diff, mindif from speed of light

# Dirty little hack to make the list of lists into a nice matrix
matcharray = ones(Int, length(matches),2)
for i in 1:length(matches)
    matcharray[i,1] = matches[i][1]
    matcharray[i,2] = matches[i][2]
end


# Sorting out which is real matches and which are fake
reals = first[matcharray[:,1], 3] .== second[matcharray[:,2], 3]
realelasticprotons = (first[matcharray[:,1], 3] .== second[matcharray[:,2], 3]) .& (second[matcharray[:,2], 2] .== "proton") .& (second[matcharray[:,2], 4] .== "hadElastic")
fakes = first[matcharray[:,1], 3] .!= second[matcharray[:,2], 3]

fig1 = scatter(first[matcharray[reals,1], 1], second[matcharray[reals,2], 1], alpha=0.1, color=:blue, label="Real hits")
scatter!(first[matcharray[fakes,1], 1], second[matcharray[fakes,2], 1], alpha=0.1, color=:red, label="Fake hits")
title!("Energy depositions of matched hits, 100MeV neutrons")
xlabel!("Energy deposition in first detector (MeV)")
ylabel!("Energy deposition in second detector (MeV)")
savefig("MatchedEnergyDepositions.svg")
display(fig1)

# Calculating the kinetic energy (assuming proton mass)
realtimes = second[matcharray[reals,2], :Column10] - first[matcharray[reals,1], :Column10]
realprotontimes = second[matcharray[realelasticprotons,2], :Column10] - first[matcharray[realelasticprotons,1], :Column10]
faketimes = second[matcharray[fakes,2], :Column10] - first[matcharray[fakes,1], :Column10]

# Here I make a cut on particles that would seem to be superluminal, for physical reasons, but in principle time uncertainty could make a particle at almost c look like it was slightly above c. I do it because sqrt of negative numbers is crashing my program
distance = 0.5025 - 0.005 # In meters. (The minus 5mm is because the protons that hit the second generate signal immediately, while the protons from the first are created uniformly through the volume by neutron scattering.)
reallorentz = 1. .- ((distance ./ realtimes)/299792458).^2
realvelocigies = (1 ./ sqrt.(reallorentz[reallorentz .>0]) .-1) * 938.27
realprotonlorentz = 1. .- ((distance ./ realprotontimes)/299792458).^2
realprotonvelocigies = (1 ./ sqrt.(realprotonlorentz[realprotonlorentz .>0]) .-1) * 938.27
fakelorentz = 1. .- ((distance ./ faketimes)/299792458).^2
fakevelocigies = (1 ./ sqrt.(fakelorentz[fakelorentz .>0]) .-1) * 938.27



bins = range(0, 250, length=250)
fig3 = stephist(realvelocigies, bins=bins, color=:black, label="Real counts")
stephist!(realprotonvelocigies, bins=bins, color=:blue, label="Real protons")
stephist!(fakevelocigies, bins=bins, color=:red, label="Fake counts")
xlabel!("Guessed incident energy (MeV)")
ylabel!("Counts")
title!("Real and fake counts, 1ns / event")
savefig("RealAndFake.svg")
display(fig3)


fig4 = scatter(second[matcharray[reals,2], 1], realvelocigies, color=:blue, alpha=0.1, markersize=1, label="Real hits")
scatter!(second[matcharray[fakes,2], 1], fakevelocigies, color=:red, alpha=0.1, markersize=1, label="Fake hits")
ylims!(0, 250)
xlabel!("Energy deposition in second detector (MeV)")
ylabel!("Guessed incident energy (MeV)")
title!("Energy deposition and guessed energies")
savefig("EnergiesAndDeposition.svg")
display(fig4)


fig5 = scatter(second[matcharray[reals,2], 1], realvelocigies, first[matcharray[reals,1], 1], color=:blue, alpha=0.1, markersize=1, label="Real hits")
scatter!(second[matcharray[fakes,2], 1], fakevelocigies, first[matcharray[fakes,1], 1], color=:red, alpha=0.1, markersize=1, label="Fake hits")
ylims!(0, 250)
xlabel!("Energy deposition in second detector (MeV)")
ylabel!("Guessed incident energy (MeV)")
zlabel!("First detector (MeV)")
title!("Energy deposition and guessed energies")
savefig("EnergiesAndDeposition.svg")
display(fig5)

#anim = Animation()
#for i in range(0, stop = 360 / 3, step = 1)
#    p = scatter(second[matcharray[reals,2], 1], realvelocigies, first[matcharray[reals,1], 1], color=:blue, alpha=0.1, markersize=1, label="Real hits", camera=(i, i/10), dpi=200)
#    scatter!(second[matcharray[fakes,2], 1], fakevelocigies, first[matcharray[fakes,1], 1], color=:red, alpha=0.1, markersize=1, label="Fake hits")
#    ylims!(0, 250)
#    xlabel!("Second detector (MeV)")
#    ylabel!("Kinetic energy (MeV)")
#    zlabel!("First detector (MeV)")
#    title!("Energy deposition and guessed energies")
#    frame(anim, p)
#end
#gif(anim, "gr.gif", fps=24)

realseconds = second[matcharray[reals,2], :]
realfirsts = first[matcharray[reals,1], :]
fakefirsts = first[matcharray[fakes,1], :]

realseconds[!,"Column11"] = realvelocigies
#sort!(realseconds, [:Column11])

realfirsts[!,"Column11"] = realvelocigies
sort!(realfirsts, [:Column9])

bins = range(0, 250, length=250)
fig6 = stephist(realfirsts[:,1] + realvelocigies, bins=bins, color=:black, label="Real counts")
stephist!(fakefirsts[:,1] + fakevelocigies, bins=bins, color=:red, label="Fake counts")
xlabel!("Proton energy @ creation (MeV)")
ylabel!("Counts")
title!("Real and fake counts, 1ns / event")

created = realfirsts[:,1] + realvelocigies
print("Values for recoil protons: \n")
print("Mean is: ", mean(created[(created .>80) .& (created .<120)]), " ", mean((realfirsts[:,7])[(realfirsts[:,7] .> 80) .&  (realfirsts[:,7] .< 120)]), "\n")
print("STD is: ", std(created[(created .>80) .& (created .<120)]), " ", std((realfirsts[:,7])[(realfirsts[:,7] .> 80) .&  (realfirsts[:,7] .< 120)]), "\n")

fitdata = created[(created .>80) .& (created .<120)]
gaussfit = fit_mle(Normal, fitdata)
print(params(gaussfit))
plot!(collect(80:120), length(fitdata) * pdf(gaussfit, collect(80:120)), color=:blue, label="ML fit")

savefig("ProtonEnergies.svg")
display(fig6)


############## Now using the angle to calculate what the incident neutron's energy must have been
a = (created ./ 938.27) .^ 2 .* cos(atan(5/50))^2  .- cos(atan(5/50))^2 .- (created ./ 938.27) .^ 2
b = 2 .* created ./ 938.27
c = - 1 .* cos(atan(5/50))^2 .* ((created/938.27) .^ 2 .-1)

incidents = (-b .- sqrt.(b.^2 .- 4 .* a .* c)) ./ (2 .* a) .* 938.27 .- 938.27
##############



bins = range(0, 250, length=250)
fig7 = stephist(incidents, bins=bins, color=:black, label="Data")
xlabel!("Determined incident neutron energy (MeV)")
ylabel!("Counts")
title!("Incident neutron energies")
print("Values for incident neutrons: \n")
print("Mean is: ", mean(incidents[(incidents .>80) .& (incidents .<120)]), "\n")
print("STD is: ", std(incidents[(incidents .>80) .& (incidents .<120)]), "\n")

fitdata = incidents[(incidents .>80) .& (incidents .<120)]
gaussfit = fit_mle(Normal, fitdata)
print(params(gaussfit))
plot!(collect(80:120), length(fitdata) * pdf(gaussfit, collect(80:120)), color=:blue, label="ML fit")

savefig("NeutronEnergies.svg")
display(fig7)


