using Plots
using Glob
using DelimitedFiles
using Statistics


filenames = glob("*.txt","SiPM test Anton/")

fig = plot(legend=false,dpi=600,minorgrid = true, minorticks=4,minorgridstyle=:solid)
risetimes = []
fwhms = []
for file in filenames
    data = readdlm(file,skipstart=6)

    max = maximum(data[:,2])

    risestart = data[findfirst(x -> x > max /10, data[:,2]),1]
    risestop = data[findfirst(x -> x > max *9 / 10, data[:,2]),1]
    push!(risetimes, risestop-risestart)

    fwhmstart = data[findfirst(x -> x > max /2, data[:,2]),1]
    fwhmstop = data[findlast(x -> x > max /2, data[:,2]),1]
    push!(fwhms,fwhmstop-fwhmstart)

    fig = plot!(data[:,1],data[:,2],alpha = 0.9)
end

print("Risetime = ", round(mean(risetimes),sigdigits=3),"ns ± ", round(std(risetimes),sigdigits=3) / sqrt(length(risetimes)),"ns \n")
print("FWHM = ",round(mean(fwhms),sigdigits=3),"ns ± ", round(std(fwhms) / sqrt(length(fwhms)),sigdigits=3),"ns \n")

xlabel!("Time (s)")
ylabel!("Voltage (V)")
xlims!(-0.00000001,0.00000002)
display(fig)
savefig(fig, "pulsesplot.png")

