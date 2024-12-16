

using Plots
using CSV
using DataFrames
using Statistics
using Distributions


# Making a function to collect energy deposit from same event but different particles (such as electron secondaries or inelastic products)
function collector(filepath)
    # Reading file with data in it
    print("\nReading: ", filepath, "\n")
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

# Function that matches hits within some time difference
function matcher(a, b, mindif, maxdif)
    aindex = 1
    matches = []
    for bindex in 1:length(b[:, end])

        # Finding the first index where the time difference is more than the minimum
        while (b[bindex, end] - a[aindex, end] > maxdif) #&& (aindex < length(a[:,end]))
            aindex += 1
        end        

        
        upperindex = aindex 
        while (b[bindex, end] - a[upperindex, end] > mindif) #&& (upperindex < length(a[:,end]))
            push!(matches, [upperindex, bindex])
            upperindex += 1
        end

        #print(aindex, "\n")
        #print(b[bindex, end])
    end


    matcharray = ones(Int, length(matches),2)
        for i in 1:length(matches)
            matcharray[i,1] = matches[i][1]
            matcharray[i,2] = matches[i][2]
        end
    return matcharray
end


# Function that takes matches and returns the energy of the incident neutron
function findincidentenergies(matcharray, first, second, distance, angle, sizecorrection)
    # Finding the time differences
    times = second[matcharray[:,2], :Column10] - first[matcharray[:,1], :Column10]

    # Calculating velocity
    traveldistance = distance - sizecorrection # In meters. (The minus 5mm is because the protons that hit the second generate signal immediately, while the protons from the first are created uniformly through the volume by neutron scattering.)
    lorentz = 1. .- ((traveldistance ./ times)/299792458).^2 # The lorentz factor
    velocigies = (1 ./ sqrt.(lorentz[lorentz .>0]) .-1) * 938.27 # Getting the velocity of protons that hit the 2nd detector, in units of MeV

    createdprotons = (first[matcharray[:,1],1])[lorentz .>0] .+ velocigies # Adding the energy deposited in the first detector to find the energy that the recoil protons are created at
    firstenergies = (first[matcharray[:,1],1])[lorentz .>0]
    secondenergies = (second[matcharray[:,2],1])[lorentz .>0]

    ############## Now using the angle to calculate what the incident neutron's energy must have been
    a = (createdprotons ./ 938.27) .^ 2 .* cos(angle)^2  .- cos(angle)^2 .- (createdprotons ./ 938.27) .^ 2
    b = 2 .* createdprotons ./ 938.27
    c = - 1 .* cos(angle)^2 .* ((createdprotons/938.27) .^ 2 .-1)


    ## Oooops, let's just do a short trick to filter those out where the discriminant D is negative, because sqrt
    anew = a[(b.^2 .- 4 .* a .*c) .>= 0]
    bnew = b[(b.^2 .- 4 .* a .*c) .>= 0]
    cnew = c[(b.^2 .- 4 .* a .*c) .>= 0]

    firstnew = firstenergies[(b.^2 .- 4 .* a .*c) .>= 0]
    secondnew = secondenergies[(b.^2 .- 4 .* a .*c) .>= 0]

    a = anew
    b = bnew
    c = cnew

    # Now solving the second order polynomial from the kinematics equation
    incidents = (-b .- sqrt.(b.^2 .- 4 .* a .* c)) ./ (2 .* a) .* 938.27 .- 938.27 

    # Returns an array with the calculated energies, as well as arrays of the deposited energies in the two detectors for the event
    return incidents, firstnew, secondnew
end



# Reading geometry file
geometry = CSV.read("geometry.txt", DataFrame; delim=",", ignorerepeated=true)

front = geometry[geometry[:,end],:]
back = geometry[.! geometry[:,end],:]

# Making a dataframe to store results
incidentframe = DataFrame(FrontIndex = Vector{Int}(), BackIndex = Vector{Int}(), Incidents = Vector{Vector{Float64}}(), FrontEnergy = Vector{Vector{Float64}}(), BackEnergy = Vector{Vector{Float64}}())

# Adding time, some nanoseconds between each beam neutron
addedtime = 20e-9

for i in 1:length(front[:,end]) # Looping through all detectors in front. Multithreaded
    # Collecting the results into detector hits instead of separate particles
    first = collector("output/front"*string(i)*".phsp")
    first[!,"Column10"] = first[:,3]*addedtime + first[:,end]
    sort!(first, [:Column10])
    for j in 1:length(back[:,end]) # Looping through all detectors in back. Multithreaded
        ############## Kinetmatics calculations
        displacementvector = [front[i,1] - back[j,1], front[i,2] - back[j,2] ,front[i,3] - back[j,3]]
        distance = sqrt(sum(displacementvector.^2))
        sizecorrection = back[j,4]/4 # The proton produces signal immediately when it hits the back scintillator
        angle = acos(sum(displacementvector .* [0,0,1]) / distance) # Calculating the angle respective to the z direction
        minimumtime = distance / 299792458 # Calculating the time it would take light to travel between the scintillators
        
        minimumenergy = 10 # Choosing a minimum energy of the travelling protons (so the lowest energy the method can see is "minimum energy + deposited energy")
        maximumtime = distance / (299792458 * sqrt(1 - (1 / (minimumenergy/938.27 +1))^2)) # Calculating the maximum time from the minimum energy
        
        ############## Reading the data files and making
        # Collecting the results into detector hits instead of separate particles

        second = collector("output/back"*string(j)*".phsp")
        second[!,"Column10"] = second[:,3]*addedtime + second[:,end]
        sort!(second, [:Column10])

        ############## Doing the time matching and energy calculation
        matchings = matcher(first, second, minimumtime, maximumtime)
        incidents, firsts, seconds = findincidentenergies(matchings, first, second, distance, angle, sizecorrection)
        push!(incidentframe, [i,j, incidents, firsts, seconds]) # Adding the energies into a dataframe
        
    end
end

print(incidentframe)

###############
# Collecting results from all individual pairs into one big list
incidents = zeros(0)
firsts = zeros(0)
seconds = zeros(0)
for i in 1:length(incidentframe[:,3])
    append!(incidents, incidentframe[i,3])
    append!(firsts, incidentframe[i,4])
    append!(seconds, incidentframe[i,5])
end


# Doing a slight bit of statistics
medi = median(incidents) # Getting a decent robust guess for the mean of the peak so I can make a un-bad cut when fitting
lower = medi - quantile(incidents, (1-0.68)/2) # Robust guesses for the standard deviation of the peak
upper = quantile(incidents, 1-(1-0.68)/2)- medi # Robust guesses for the standard deviation of the peak

# ML fit, cutting data 20% below and above the calculated mean
fitdata = incidents[(incidents .> medi - lower*3) .& (incidents .< medi + upper*3)] # Cutting a roughly 3sigma region around the peak
gaussfit = fit_mle(Normal, fitdata)
print(params(gaussfit), (upper + lower)/2, "\n")

fig1 = plot()
for i in 1:length(incidentframe[:,3])
    stephist!(incidentframe[i,3], label="Front"*string(incidentframe[i,1])*", Back"*string(incidentframe[i,2]))
    #print(length(incidentframe[i,3]),"\n")
end

xlims!(0,medi*2)
savefig("output/plots/SeparateEnergies.svg")
title!("Energies from detector pairs")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("NeutronHeatmap.svg")
display(fig1)

fig2 = histogram(incidents[incidents .< medi*2])
title!("Total energy spectrum")
xlabel!("Energy (MeV)")
ylabel!("Counts")
savefig("NeutronHeatmap.svg")
savefig("plots/TotalEnergies.svg")
display(fig2)


fig3 = histogram2d(incidents[incidents .< medi*2], firsts[incidents .< medi*2])
title!("Incident energy and first detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in first detector (MeV)")
savefig("plots/FirstHeatmap.svg")
display(fig3)

fig4 = histogram2d(incidents[incidents .< medi*2], seconds[incidents .< medi*2])
title!("Incident energy and second detector")
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Energy in second detector (MeV)")
savefig("plots/SecondHeatmap.svg")
display(fig4)

# Writing into results file
file = open("output/results.csv", "a")
write(file, string(params(gaussfit)[1]) * ", " * string(params(gaussfit)[2]))
close(file)

