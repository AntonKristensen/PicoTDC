
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Glob

include("../cutting.jl")

open("stochasticresults.csv", "w") do f
    write(f, "addedtime (ps), correct, fake, S/N >25MeV")
end

files = glob("output/matches*")
for filepath in files
    firstread = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)
    
    lines = readlines(filepath)[2:end]
    # Remove lines if they wrongly contain header
    deleteat!(lines, firstread[:,1] .== "Incident energy")

    # Write the remaining lines back to the file
    open(filepath, "w") do f
        write(f, "Incident energy, first deposit (MeV), second deposit (MeV), front detector, back detector, front event, back event\n")
        write(f, join(lines, "\n"))
    end
    
    noncutdata = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)

    #cleaneddata.age = parse.(Int, df.age)
    data = cut(noncutdata)

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



    geometry = CSV.read("geometry.txt", DataFrame; delim=",", ignorerepeated=true, ignoreemptyrows=true)

    # Making it so that the neutron beam is wide enough to cover all the front scintillators. 
    area = 10000 *  maximum(abs.(geometry[geometry[:,end], 1]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))) * maximum(abs.(geometry[geometry[:,end], 2]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2)))

    frontmaxx = string(maximum(abs.(geometry[geometry[:,end], 1]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
    frontmaxy = string(maximum(abs.(geometry[geometry[:,end], 2]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
    frontmaxsize = string(maximum(abs.(geometry[geometry[:,end], 4]./2)))

    println(area)
    flux = 1 /( 1.0e-12 * parse(Float64, filepath[15:end-4]) * area)

    correct = length(cincidents)
    fake =  length(fincidents)

    fig1 = histogram(cincidents, bins=0:1:maximum(cincidents), color=:black, alpha=0.5, label="Correct events")
    histogram!(fincidents, bins=0:1:maximum(cincidents), color=:red, alpha=0.5, label="Stochastic events")
    title!(string(round(flux, sigdigits=3)) * " neutrons s⁻¹ cm⁻²")
    xlabel!("Calculated neutron energy (MeV)")
    ylabel!("Counts")
    display(fig1)
    savefig("plots/" * filepath[15:end-4] * ".svg")


    println(filepath, ": Correct: ", correct, ", Fake: ", fake, ", Fake above 25 MeV: ", sum(fincidents .> 25))


    open("stochasticresults.csv", "a+") do f
        write(f, "\n",filepath[15:end-4], ", ", string(correct), ", ", string(fake), ", ", string(sum(cincidents .> 25) / sum(fincidents .> 25)))
    end
end
