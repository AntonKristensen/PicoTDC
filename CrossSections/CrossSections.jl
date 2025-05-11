

# Importing packages
using Plots # Plotting
using CSV # Allows CSV.read(filepath, DataFrame; header=0, delim=",")
using DataFrames # DataFrame type for use with CSV
using Statistics # Statistics stuff
using StatsBase # Can compute histograms by fit(Histogram, data, bins)
using LinearAlgebra # Has the normalize function for histograms
using Distributions # For doing statistics with like gaussians etc
using LegendrePolynomials # Computing Pl(cos, l)
using Interpolations # Easy linear interpolation
using NumericalIntegration # For integrating: integrate(x,y)

# Loading data file
function crosssectionloader(filepath)
    data = CSV.read(filepath, DataFrame; header=0, delim=" ", ignorerepeated=true)
    data = data[:, 1:6] # Throwing away metadata
    data = replace.(data, "-" => "e-") # Fixing scientific notation
    data = replace.(data, "+" => "e+") # Same
    data = parse.(Float64, data) # Turn it from strings to floats

    # Folding out from FORTRAN style data tape structure
    energies = zeros(length(data[:,1]) * 3)
    sections = zeros(length(data[:,1]) * 3)
    for n in 1:length(data[:,1])
        energies[n*3 - 2] = data[n,1]
        energies[n*3 - 1] = data[n,3]
        energies[n*3] = data[n,5]

        sections[n*3 - 2] = data[n,2]
        sections[n*3 - 1] = data[n,4]
        sections[n*3] = data[n,6]
    end
    return energies ./ 1e6, sections
end

totale, totals = crosssectionloader("total.txt")
elastice, elastics = crosssectionloader("elastic.txt")

fig1 = plot(totale, totals, xaxis=:log, yaxis=:log, label = "Total", color=:blue, linewidth=3, alpha=0.5, grid=true, minorgrid=true)
plot!(elastice, elastics, label="Elastic", color=:red, linewidth=3, alpha=0.5)
title!("Neutron-proton cross sections")
xlims!(1, 250)
xaxis!(minorticks=10)
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Cross section (barn)")
savefig("total.png")
display(fig1)

function adde(string) # OMG I hate this dataformat. I have to write this function to check where the exponential starts
    if isdigit(string[end-1]) # Checks if there is one or two numbers in the exponential
        result = string[1:end-3] * "e" * string[end-2:end]
    else
        result = string[1:end-2] * "e" * string[end-1:end]
    end
    return result
end

function angularloader(filepath)
    file = open(filepath, "r") # Opening the data file
    energies = zeros(0)
    coeffs = DataFrame([0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.], :auto) # empty dataframe to put coefficients into
    while !eof(file) # Reads until the file ends
        line = readline(file)
        if line[1:11] == " 0.000000+0"
            push!(energies, parse(Float64, adde(line[12:22]))) # Add the energy

            number = trunc(Int, parse(Int, line[54:55])/6)

            newline = readline(file)

            loworder = zeros(6)

            # Reading the individual legendre coeffs
            #loworder[1] = parse(Float64, newline[1:8] * "e" * newline[9:11])
            loworder[1] = parse(Float64, adde(newline[1:11]))
            loworder[2] = parse(Float64, adde(newline[12:22]))
            loworder[3] = parse(Float64, adde(newline[23:33]))
            loworder[4] = parse(Float64, adde(newline[34:44]))
            loworder[5] = parse(Float64, adde(newline[45:55]))
            loworder[6] = parse(Float64, adde(newline[56:66]))

            if number == 2
                newerline = readline(file)
                highorder = zeros(6)
                highorder[1] = parse(Float64, adde(newerline[1:11]))
                highorder[2] = parse(Float64, adde(newerline[12:22]))
                highorder[3] = parse(Float64, adde(newerline[23:33]))
                highorder[4] = parse(Float64, adde(newerline[34:44]))
                highorder[5] = parse(Float64, adde(newerline[45:55]))
                highorder[6] = parse(Float64, adde(newerline[56:66]))
            else
                highorder = zeros(6)
            end
            result = [loworder ; highorder]
            push!(coeffs, result)
        end
        #break
    end
    close(file) # Making sure to close the file
    deleteat!(coeffs, 1)
    return energies / 1.0e6, coeffs
end

e, c = angularloader("angular.txt")

fig2 = plot(e, abs.(c[:,1]), label="l=1", xaxis=:log, yaxis=:log, linewidth=3, alpha=0.5, grid=true, minorgrid=true, legend=:topleft)
for n in 2:length(c[1,:])
    plot!(e[abs.(c[:,n]) .> 0], abs.((c[:,n])[abs.(c[:,n]) .> 0]), label="l="*string(n), linewidth=3, alpha=0.5)
end
#plot!(e, c[:,2], label="l=1", linewidth=3, alpha=0.5)
title!("Neutron-proton legendre coefficients")
#xlims!(minimum(e), 250)
xlims!(1, 250)
ylims!(1e-7, 1)
xaxis!(minorticks=10)
xlabel!("Energy of incident neutron (MeV)")
ylabel!("Legendre polynomial coefficient")
savefig("Legendre.png")
display(fig2)


#############################################

# Function for finding the value of the probability function at some angle and some index of the energy list
function evaluateatcos(cos, energyindex)
    lfactor = [(2*l +1)/2 for l in 1:12]
    return sum((collectPl(cos, lmax=12))[1:end] .* collect(c[energyindex,:]) .* lfactor) + 0.5
end

cosses = collect(range(-1, 1, length=500))
probs = evaluateatcos.(cosses, length(e))
fig3 = plot(cosses, probs, label="E=200MeV", linewidth=3, alpha=1, grid=true, minorgrid=true, legend=:topleft)
plot!(cosses, evaluateatcos.(cosses, findfirst(x -> x >= 100, e)), label="E="*string(e[findfirst(x -> x > 100, e)])*"MeV", linewidth=3, alpha=1)
plot!(cosses, evaluateatcos.(cosses, findfirst(x -> x >= 60, e)), label="E="*string(e[findfirst(x -> x > 60, e)])*"MeV", linewidth=3, alpha=1)
plot!(cosses, evaluateatcos.(cosses, findfirst(x -> x >= 21, e)), label="E="*string(e[findfirst(x -> x > 25, e)])*"MeV", linewidth=3, alpha=1)
#plot!(cosses, evaluateatcos.(cosses, findfirst(x -> x >= 10, e)), label="E="*string(e[findfirst(x -> x > 10, e)])*"MeV", linewidth=3, alpha=0.5)
#plot!(cosses, evaluateatcos.(cosses, findfirst(x -> x >= 1, e)), label="E="*string(e[findfirst(x -> x > 1, e)])*"MeV", linewidth=3, alpha=0.5)
title!("JENDL angular distribution of n-p scattering")
xlims!(0,1)
ylims!(0,1.3)
xlabel!("Lab frame cos(θ)")
ylabel!("Relative scattering probability")
savefig("angular.pdf")
display(fig3)

#################


function argonneloader(filepath)
    data1 = CSV.read(filepath, DataFrame; header=5, footerskip=15, delim=" ", ignorerepeated=true)
    data2 = CSV.read(filepath, DataFrame; header=20, footerskip=0, delim=" ", ignorerepeated=true)
    
    return innerjoin(data1, data2, on = :Tlab)
end
argonneshifts = argonneloader("argonne.txt")
#print(argonneshifts)



function nijmloader(filepath)
    data = CSV.read(filepath, DataFrame; header=5, delim=" ", ignorerepeated=true)
    
    return data
end
nijmpartial = nijmloader("nijm.txt")[1:4:end,:]
nijmnonlocal = nijmloader("nijm.txt")[2:4:end,:]
nijmlocal = nijmloader("nijm.txt")[3:4:end,:]
reid = nijmloader("nijm.txt")[4:4:end,:]

#print(nijmpartial)

# Using the Argonne 18 potential:





function CoMtoLab(angle, energy)
    pmass = 938.272 # In MeV/c²
    pl = sqrt((pmass + energy)^2 - pmass^2) # momentum of projectile particle in lab frame
    v = pl / (energy + pmass) # velocity of lorentz transformation
    p = (pl - v * energy)/sqrt(1-v^2) # momentum in CoM
    u = p / sqrt(pmass^2 + p^2) # Velocity of projectile in CoM

    eprime = pmass /sqrt(1 - u^2)
    qprime = pmass * u /sqrt(1-u^2)

    radians = angle#acos(cosangle)

    tan = sqrt(1-v^2) * qprime * sin(radians) / (qprime*cos(radians) + v * eprime) # tan of angle in lab

    return [tan, v]
end


function argonne(angle, energy, phaseshifts=[7.973, -1.368, -20.793, -21.257, 7.323, 24.194], format="") # Function that returns the scattering amplitude of the argonne potential for n-p scattering. Units for energy are MeV
    # cosangle is the cosine of the angle in the CoM frame, and energy is the energy in the lab frame
    # The three different l contributions

    phaseshifts = Vector(phaseshifts)
    if format == "argonne"
        l0 = phaseshifts[[2,12]]
        l1 = phaseshifts[[4,5,6,10]]
        l2 = phaseshifts[[3,14,15,16]]
        l3 = phaseshifts[[8,9,11]]
        l4 = phaseshifts[[18]]
    elseif format == "nijm"
        l0 = phaseshifts[[2,6]]
        l1 = phaseshifts[[3,4,5,11]]
        l2 = phaseshifts[[8,9,10]]
        l3 = phaseshifts[[13]]
        l4 = phaseshifts[[]]
    else
        #test = [7.973, -1.368, -20.793, -21.257, 7.323, 24.194]
        l0 = phaseshifts[1:3]
        l1 = phaseshifts[4:5]
        l2 = phaseshifts[5:6]
        l3 = [0]
        l4 = [0]
    end

    cosangle = cos(angle)
    zero = sum(  (2*0 +1)  *  exp.(im * l0)  .* sin.(2*pi * l0 / 360)  *  Pl(cosangle, 0)  )
    one = sum(  (2*1 +1)  *  exp.(im * l1)  .* sin.(2*pi * l1 / 360)  *  Pl(cosangle, 1)  )
    two = sum(  (2*2 +1)  *  exp.(im * l2)  .* sin.(2*pi * l2 / 360)  *  Pl(cosangle, 2)  )
    three = sum(  (2*3 +1)  *  exp.(im * l3)  .* sin.(2*pi * l3 / 360)  *  Pl(cosangle, 3)  )
    four = sum(  (2*4 +1)  *  exp.(im * l4)  .* sin.(2*pi * l4 / 360)  *  Pl(cosangle, 4)  )

    hbar = 5.582e-22 # MeV s
    reducedmass = 469.459 # MeV
    c = 299792458 # m/s

    comenergy = energy * CoMtoLab(cosangle, energy)[2] # Use the energy in CoM frame, not lab frame
    k = sqrt(2 * reducedmass/(c^2) * comenergy / (hbar^2)) # units of m⁻¹

    result = abs2((zero + one + two + three + four) / k) # Norm squaring

    tanangle = CoMtoLab(angle, energy)[1]
    # Time for trigonemtry to find which quadrant the angle is in
    if tanangle >= 0 && cosangle >= 0
        labangle = atan(tanangle)
    elseif tanangle < 0 && cosangle < 0
        labangle = pi +atan(tanangle)
    elseif tanangle >= 0 && cosangle < 0
        labangle = atan(tanangle)
    elseif tanangle < 0 && cosangle >= 0
        labangle = pi +atan(tanangle)
    end

    return result * 1e28, labangle # Units of barn, also returns the tangent of the lab angle
end

angles = collect(range(0,  2*pi, length=500))
argonnes = zeros(length(angles))
nijms = zeros(length(angles))
nonlocals = zeros(length(angles))
locals = zeros(length(angles))
reids = zeros(length(angles))
labangles = zeros(length(angles))
for n in 1:length(angles)
    argonnes[n], labangles[n] = argonne(angles[n],argonneshifts[8,1], argonneshifts[8,:], "argonne")
    nijms[n], labangles[n] = argonne(angles[n],nijmpartial[8,1], nijmpartial[8,:], "nijm")
    nonlocals[n], labangles[n] = argonne(angles[n],nijmnonlocal[8,1], nijmnonlocal[8,:], "nijm")
    locals[n], labangles[n] = argonne(angles[n],nijmlocal[8,1], nijmlocal[8,:], "nijm")
    reids[n], labangles[n] = argonne(angles[n],reid[8,1], reid[8,:], "nijm")



end


fig4 = plot(cos.(labangles), argonnes, label="200MeV Argonne", linewidth=3, alpha=1, grid=true, minorgrid=true, legend=:topleft)
plot!(cos.(labangles), nijms, label="200MeV Nijm", linewidth=3, alpha=1)
plot!(cos.(labangles), nonlocals, label="200MeV Nijm non-local", linewidth=3, alpha=1)
plot!(cos.(labangles), locals, label="200MeV Nijm local", linewidth=3, alpha=1)
plot!(cos.(labangles), reids, label="200MeV Reid", linewidth=3, alpha=1)

title!("Differential cross section from Argonne potential")
xlims!(0,1)
ylabel!("dσ/dΩ (barns)")
xlabel!("Lab frame cos(θ)")
savefig("argonne.pdf")
display(fig4)


#################### From topas simulation
# Reading data
topasdata = CSV.read("simulated/100mil.phsp", DataFrame; header=false, delim=" ", ignorerepeated=true)

x = topasdata[:,1]
y = topasdata[:,2]
z = topasdata[:,3]

phis = atan.(y ./ x)
thetas = atan.(sqrt.(y.^2 .+ x.^2) ./ z)

angularhist = histogram(cos.(thetas), bins=200, label="200MeV n-p")
title!("Geant4 n-p scattering")
xlims!(0,1)
xlabel!("Lab frame cos(θ)")
ylabel!("Counts")
savefig("simulation.pdf")
display(angularhist)


############## Making a combined plot, so needing some normalisation
intervalues = 0:1/100:1

an = cos.(labangles)[cos.(labangles) .>= 0]
nij = nijms[cos.(labangles) .>= 0]
nij = nij ./ abs(integrate(an,nij))

jcosses = cosses[cosses .>= 0]
jprobs = probs[cosses .>= 0]
jprobs = jprobs ./ abs(integrate(jcosses,jprobs))

ghist = fit(Histogram, cos.(thetas), collect(0:1/100:1))
ghist = normalize(ghist, mode=:pdf)

#combinedplot = stephist(cos.(thetas), bins=100, linewidth=3, normalize=:pdf, label="Geant4")
combinedplot = plot(collect(0+1/100:1/100:1), ghist.weights, linewidth=3, label="Geant4", grid=true, minorgrid=true)
plot!(an, nij, label="Nijm", linewidth=3, alpha=1)
plot!(jcosses, jprobs, label="JENDL", linewidth=3)
title!("n-p scattering at 200MeV")
xlabel!("Lab frame cos(θ)")
ylabel!("Relative angular cross section")
savefig("Combined.pdf")
display(combinedplot)
