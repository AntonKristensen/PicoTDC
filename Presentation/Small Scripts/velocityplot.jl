
using Plots


x = collect(0:0.1:500) # List of energies

function velocity(E; mass=938.272) # Function for finding velocity of proton/neutron from kinetic energy. Using E in units of MeV, v as ratio to c
    return sqrt.(1 .- (1 ./ (E / mass .+1)).^2)
end

fig = plot(x, velocity(x), linewidth=5, label="Proton")
plot!(x, velocity(x, mass=80*938.272), linewidth=5, label="A=80")
title!("Velocity of proton as function of kinetic energy")
xlabel!("Kinetic energy (MeV)")
ylabel!("Velocity (c)")
savefig("Velocityplot.png")
display(fig)


println("100MeV: ", velocity(100, mass=80*938.272),"c. 250MeV: ", velocity(250, mass=80*938.272), "c.")

timedif = 1 / (velocity(100, mass=80*938.272) * 299792458) - 1 / (velocity(200, mass=80*938.272) * 299792458)
println(timedif)