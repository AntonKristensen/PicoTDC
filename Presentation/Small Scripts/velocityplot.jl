
using Plots


x = collect(0:0.1:500) # List of energies

function velocity(E) # Function for finding velocity of proton/neutron from kinetic energy. Using E in units of MeV, v as ratio to c
    return sqrt.(1 .- (1 ./ (E / 938.272 .+1)).^2)
end

fig = plot(x, velocity(x), linewidth=5)
title!("Velocity of proton as function of kinetic energy")
xlabel!("Kinetic energy (MeV)")
ylabel!("Velocity (c)")
savefig("Velocityplot.png")
display(fig)


println("20MeV: ", velocity(20),"c. 250MeV: ", velocity(250), "c.")

timedif = 1 / (velocity(20) * 299792458) - 1 / (velocity(250) * 299792458)
println(timedif)