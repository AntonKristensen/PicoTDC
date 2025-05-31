
using Plots
using CSV
using DataFrames
using Statistics
using Distributions


function reader(filepath)
    return CSV.read(filepath, DataFrame; header=false, delim=" ", ignorerepeated=true)
end

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
function matcher(a, b, mindif, maxdif) #Has time cutoff energy as keyword argument, with default being 0
    aindex = 1
    matches = []
    for bindex in 1:length(b[:, end])

        # Finding the first index where the time difference is more than the minimum
	    while (b[bindex, end] - a[aindex, end] > maxdif) && (aindex < length(a[:,end]) -1 ) 
            aindex += 1
        end        

        
        upperindex = aindex 
        while (b[bindex, end] - a[upperindex, end] > mindif) && (upperindex < length(a[:,end]) -1)
            push!(matches, [upperindex, bindex])
            upperindex += 1
        end

    end

    matcharray = ones(Int, length(matches),2)
        for i in 1:length(matches)
            matcharray[i,1] = matches[i][1]
            matcharray[i,2] = matches[i][2]
        end
    return matcharray
end


# Function that takes matches and returns the energy of the incident neutron
function findincidentenergies(matcharray, first, second, distance, angle, sizecorrection; timeuncertainty=0)
    # Finding the time differences
    times = (second[matcharray[:,2], :Column10] .+ timeuncertainty) .- (first[matcharray[:,1], :Column10] .+ timeuncertainty)
    #times = second[matcharray[:,2], :Column10] - first[matcharray[:,1], :Column10]


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
    incidenthadrons = (-b .- sqrt.(b.^2 .- 4 .* a .* c)) ./ (2 .* a) .* 938.27 .- 938.27 

    # Returns an array with the calculated energies, as well as arrays of the deposited energies in the two detectors for the event
    return incidenthadrons, firstnew, secondnew
end


function seedrng(seed) # Function for making a random number, but using a specific seed, so I can use G4event number as seed, to make sure that the same amount of random time is added to events in different detectors
    Random.seed!(seed) # Manually sets the random seed to the input, ensuring determinism in the Reading
    return randn() # Returns a normally distributed random number
end


function detectorlooping(; geofile = "geometry.txt", addedtime = 1.e-9, minimumenergy = 1., timeuncertainty = 0., threshold=0) # Function for looping through all the detector outputs and matching events 
    # Reading geometry file
    geometry = CSV.read(geofile, DataFrame; delim=",", ignorerepeated=true)

    front = geometry[geometry[:,end],:]
    back = geometry[.! geometry[:,end],:]  

    # Making a dataframe to store results
    incidentframe = DataFrame(FrontIndex = Vector{Int}(), BackIndex = Vector{Int}(), Incidents = Vector{Vector{Float64}}(), FrontEnergy = Vector{Vector{Float64}}(), BackEnergy = Vector{Vector{Float64}}(), FrontG4Event = FrontEnergy = Vector{Vector{Int}}(), BackG4Event = FrontEnergy = Vector{Vector{Int}}())
    numberofpairs = length(front[:,end]) * length(back[:,end])
    defaultframe = copy(incidentframe)
    push!(defaultframe, [0, 0, [], [], [], [], []])
    for n in 1:numberofpairs
        append!(incidentframe, defaultframe)
    end

    Threads.@threads for i in 1:length(front[:,end]) # Looping through all detectors in front. Multithreaded
        if filesize("output/front"*string(i)*".phsp") != 0 # Don't do it if the file is empty
            # Collecting the results into detector hits instead of separate particles
            first = collector("output/front"*string(i)*".phsp")
            #first[!,"Column10"] = (first[:,3]*addedtime + first[:,end] ) + randn(length(first[:,end])) * timeuncertainty # Adding time, some nanoseconds between each beam neutron. Consider making this more sophisticated, so that the times are distributed randomly according to some distribution.
            first[!,"Column10"] = (first[:,3]*addedtime + first[:,end] ) + randn(length(first[:,end])) * timeuncertainty # Adding time, some nanoseconds between each beam neutron. Consider making this more sophisticated, so that the times are distributed randomly according to some distribution.
            
            for n in 1:length(first[:,3])
                first[n,"Column10"] = first[n,"Column10"] + addedtime * seedrng(first[n, 3]) # Add a random amount of time, but needs to be same for same G4event in different detectors
            end
            
            sort!(first, [:Column10])
            #first = first[first[:,1] .> threshold ,:] # Sets a lower energy deposition limit
            Threads.@threads for j in 1:length(back[:,end]) # Looping through all detectors in back. Multithreaded
                if filesize("output/back"*string(j)*".phsp") != 0 # Don't do it if the file is empty   
                
                    ############## Kinetmatics calculations
                    displacementvector = [front[i,1] - back[j,1], front[i,2] - back[j,2] ,front[i,3] - back[j,3]]
                    distance = sqrt(sum(displacementvector.^2))
                    sizecorrection = back[j,4]/4 # The proton produces signal immediately when it hits the back scintillator
                    angle = acos(sum(displacementvector .* [0,0,1]) / distance) # Calculating the angle respective to the z direction
                    minimumtime = distance / 299792458 # Calculating the time it would take light to travel between the scintillators
                    
                    # Choosing a minimum energy of the travelling protons (so the lowest energy the method can see is "minimum energy + deposited energy")
                    maximumtime = distance / (299792458 * sqrt(1 - (1 / (minimumenergy/938.27 +1))^2)) # Calculating the maximum time from the minimum energy
                    
                    ############## Reading the data files and making
                    # Collecting the results into detector hits instead of separate particles
                    second = collector("output/back"*string(j)*".phsp")
                    #second[!,"Column10"] = second[:,3]*addedtime + second[:,end] + randn(length(second[:,end])) * timeuncertainty  # Adding some time to each event. Consider making this more sophisticated, so that the times are distributed randomly according to some distribution.
                    second[!,"Column10"] = second[:,3]*addedtime + second[:,end] + randn(length(second[:,end])) * timeuncertainty  # Adding some time to each event. Consider making this more sophisticated, so that the times are distributed randomly according to some distribution.
                    for n in 1:length(first[:,3])
                        second[n,"Column10"] = second[n,"Column10"] + addedtime * seedrng(second[n, 3]) # Add a random amount of time, but needs to be same for same G4event in different detectors
                    end
                    #second = second[second[:,1] .< threshold ,:] # Sets a lower energy deposition limit
                    sort!(second, [:Column10])
                    ############## Doing the time matching and energy calculation
                    matchings = matcher(first, second, minimumtime, maximumtime)
                    incidents, firsts, seconds = findincidentenergies(matchings, first, second, distance, angle, sizecorrection)
                    firstevents = first[matchings[:,1], 3]
		    secondevents = second[matchings[:,2],3]

                    #push!(incidentframe, [i,j, incidents, firsts, seconds]) # Adding the energies into a dataframe
                    incidentframe[(i-1)*length(front[:,end]) + j, :] = [i,j, incidents, firsts, seconds, firstevents, secondevents]
                end
            end
        end
    end
    return incidentframe
end



# Writing the matches into a file
function matchwriter(dataframe; file = "output/matches.csv")
    matchfile = open(file, "a")
    if filesize("output/matches.csv") == 0 # Checks if the file is empty, and writes a header if it is
        write(matchfile, "Incident energy, first deposit (MeV), second deposit (MeV), front detector, back detector, front event, back event \n")
    end
    for i in 1:length(dataframe[:,1])
        if length(dataframe[i,3]) > 0 # Don't write if there ain't no data
            for j in 1:length(dataframe[i,3])
                write(matchfile, string(dataframe[i,3][j]) * ", " * string(dataframe[i,4][j]) * ", " * string(dataframe[i,5][j]) * ", " * string(dataframe[i,1]) *  ", " * string(dataframe[i,2]) * "," *  string(dataframe[i,6][j]) * ", " * string(dataframe[i,7][j]) *   "\n")   
            end
        end
    end
    flush(matchfile) # Ensures that the data is written before moving on, to avoid writing too many headers for example
    close(matchfile)
end
