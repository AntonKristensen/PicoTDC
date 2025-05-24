using Plots

function energycalculation(energy, angle)
    # Coefficients of the polynomial equation
    a = (energy ./ 938.27) .^ 2 .* cos(angle)^2  .- cos(angle)^2 .- (energy ./ 938.27) .^ 2
    b = 2 .* energy ./ 938.27
    c = - 1 .* cos(angle)^2 .* ((energy/938.27) .^ 2 .-1)

    # Now solving the second order polynomial from the kinematics equation
    resultenergy = (-b .- sqrt.(b.^2 .- 4 .* a .* c)) ./ (2 .* a) .* 938.27 .- 938.27 

    return resultenergy
end

function rad(degree)
    return degree * 2 * pi / 360
end

function deg(radian)
    return radian * 360 / (2 * pi)
end

function gibenergies(energy, centerangle)
    anglecollection = centerangle-2:0.2:centerangle+2
    anglers = rad(collect(anglecollection))
    energies = zeros(length(anglers))
    for n in 1:length(anglers)
        energies[n] = energycalculation(energy, anglers[n])
    end
    return deg(anglers) .- centerangle, (energies .- energycalculation(energy, rad(centerangle))) ./ energycalculation(energy, rad(centerangle)) * 100
end



fig1 = plot(gibenergies(50, 2), color=:black, label="50 MeV, 2 degrees")
plot!(gibenergies(50, 5), color=:blue, label="50 MeV, 5 degrees")
plot!(gibenergies(50, 15), color=:red, label="50 MeV, 15 degrees")
plot!(gibenergies(200, 2), color=:black, linestyle=:dash, label="200 MeV, 2 degrees")
plot!(gibenergies(200, 5), color=:blue, linestyle=:dash, label="200 MeV, 5 degrees")
plot!(gibenergies(200, 15), color=:red, linestyle=:dash, label="200 MeV, 15 degrees")
title!("Energy error from angle error")
xlabel!("Error in angle (degrees)")
ylabel!("Relative error in energy (%)")
savefig("AngleError.svg")
display(fig1)
