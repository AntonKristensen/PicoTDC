
using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Glob
###########
function meaner(datalist)
    results = []
    for i in 1:length(datalist[1])
        values = []
        for j in 1:length(datalist)
            push!(values, datalist[j][i])
        end
        push!(results, mean(values))
    end
    return results
end

function deviater(datalist)
    results = []
    for i in 1:length(datalist[1])
        values = []
        for j in 1:length(datalist)
            push!(values, datalist[j][i])
        end
        push!(results, std(values)./ sqrt(length(values))) # Calculate the uncertainty as std/sqrt(N)
    end
    return results
end
###########

medians = []
spreads = []
fiddimedians = []
fiddispreads = []
hunnimedians = []
hunnispreads = []
for filepath in glob("results*")
    data = CSV.read(filepath, DataFrame; header=1, delim=",", ignorerepeated=false)
    push!(medians, data[:,1])
    push!(spreads, data[:,2])
    push!(fiddimedians, data[:,3])
    push!(fiddispreads, data[:,4])
    push!(hunnimedians, data[:,5])
    push!(hunnispreads, data[:,6])
end


uncertaintyplot = plot(meaner(medians), 100 * meaner(spreads) ./meaner(medians), yerr=100 * deviater(spreads)./meaner(medians), color=:black, markerstrokecolor=:black,label="Ideal", size=(600,400))
plot!(meaner(fiddimedians), 100 * meaner(fiddispreads)./meaner(fiddimedians), yerr=100 * deviater(fiddispreads)./meaner(fiddimedians), color=:blue, markerstrokecolor=:blue,label="50ps")
plot!(meaner(hunnimedians), 100 * meaner(hunnispreads)./meaner(hunnimedians), yerr=100 * deviater(hunnispreads)./meaner(hunnimedians), color=:red, markerstrokecolor=:red, label="100ps")
title!("Energy Uncertainties")
xlabel!("Neutron energy (MeV)")
ylabel!("Uncertainty (%)")
savefig("plots/uncertainties.png")
display(uncertaintyplot)

