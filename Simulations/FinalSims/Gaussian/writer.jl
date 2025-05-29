
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Random
################## Writing overall Geant4 control
file = open("output/neutrons.topas", "a")
rng = string(abs(rand(Int32)))
write(file, "i:Ts/Seed = " * rng * "\n")



particles = ARGS[1]
energy = ARGS[2]

write(file, "ic:So/Beam/NumberOfHistoriesInRun   = " * string(particles)* "\n")


################### Handling detector geometry
geometry = CSV.read("geometry.txt", DataFrame; delim=",", ignorerepeated=true, ignoreemptyrows=true)

# Making it so that the neutron beam is wide enough to cover all the front scintillators. 
frontmaxx = string(maximum(abs.(geometry[geometry[:,end], 1]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
frontmaxy = string(maximum(abs.(geometry[geometry[:,end], 2]) .+ maximum(abs.(geometry[geometry[:,end], 4]./2))))
frontmaxsize = string(maximum(abs.(geometry[geometry[:,end], 4]./2)))

write(file, 
"dc:So/Beam/BeamPositionCutoffX      = " *frontmaxx* " m
dc:So/Beam/BeamPositionCutoffY      = " *frontmaxy* " m\n\n"*
"dv:So/Beam/BeamEnergySpectrumValues    = 4 90 91 109 110  MeV \n"*
"uv:So/Beam/BeamEnergySpectrumWeights   = 4 0 1 1 0 \n"
)


# Making it so the world is only slightly larger enough to contain the components
maxx = string(maximum(abs.(geometry[:,1])) + maximum(abs.(geometry[:,4])) + 0.01)
maxy = string(maximum(abs.(geometry[:,2])) + maximum(abs.(geometry[:,4])) + 0.01)
maxz = string(maximum(abs.(geometry[:,3])) + maximum(abs.(geometry[:,4])) + 0.01)

write(file, "d:Ge/World/HLX = "*maxx*" m
d:Ge/World/HLY = "*maxy*" m
d:Ge/World/HLZ = "*maxz*" m
s:Ge/World/Material = \"Air\"\n\n)

#s:Ge/World/Material = \"Vacuum\"\n\n")

frontcount=0
backcount=0


for i in 1:length(geometry[:,1])
    # Counting front and back detectors for naming scheme
    if geometry[i,end]
        global frontcount += 1
        name = "front"*string(frontcount)
    else
        global backcount += 1
        name = "back"*string(backcount)
    end

    # Reading geometry parameters
    size = string(geometry[i,4]/2)
    x = string(geometry[i,1])
    y = string(geometry[i,2])
    z = string(geometry[i,3])

    write(file,"s:Ge/"*name*"/Type     = \"TsBox\"
s:Ge/"*name*"/Material = \"Plastic\"
s:Ge/"*name*"/Parent   = \"World\"
d:Ge/"*name*"/HLX      = "*size*" m
d:Ge/"*name*"/HLY      = "*size*" m
d:Ge/"*name*"/HLZ      = "*size*" m
d:Ge/"*name*"/TransX   = "*x*" m
d:Ge/"*name*"/TransY   = "*y*" m
d:Ge/"*name*"/TransZ   = "*z*" m
s:Sc/"*name*"/Quantity                  = \"rindom\"
s:Sc/"*name*"/Component                 = \""*name*"\"
s:Sc/"*name*"/OutputType 		     = \"ASCII\"
b:Sc/"*name*"/OutputToConsole           = \"TRUE\"
s:Sc/"*name*"/OutputFile = \"output/"*name*"\"
s:Sc/"*name*"/IfOutputFileAlreadyExists = \"Overwrite\"\n\n")

end

if ARGS[3] == "1"
    write(file, "s:Gr/ViewA/Type              = \"OpenGL\"
    Ts/UseQt = \"True\" #This one makes user interface
    b:Ts/PauseBeforeQuit = \"True\"\n\n")
end

close(file)

