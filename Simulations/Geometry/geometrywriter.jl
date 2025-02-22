using Plots
using CSV
using DataFrames
using Statistics
using Distributions
using Random


# Number of scintillator detectors

detectors =parse(Int, ARGS[1])
size = 0.006
d = 0.01
zdif = 0.25


number = ceil(sqrt(detectors+1))-1 # "ring" number

xs = []
ys = []
fronts = []
# Loop to make detectors
for n in 1:number
    # Start at x=n and y=n (right corner)
    x = n
    y = n
    front = true

    # Going down
    while y > -n
        push!(xs, x)
        push!(ys, y)
        push!(fronts, front)
        front = !front
        y-=1
    end

    # Going left
    while x > -n
        push!(xs, x)
        push!(ys, y)
        push!(fronts, front)
        front = !front
        x-=1
    end

    # Going up
    while y < n
        push!(xs, x)
        push!(ys, y)
        push!(fronts, front)
        front = !front
        y+=1
    end

    # Going right
    while x < n
        push!(xs, x)
        push!(ys, y)
        push!(fronts, front)
        front = !front
        x+=1
    end
end

# Gathering results in lists
fx = xs[fronts .== true]
fx = fx[1:detectors]
fy = ys[fronts .== true]
fy = fy[1:detectors]
bx = xs[fronts .== false]
bx = bx[1:detectors]
by = ys[fronts .== false]
by = by[1:detectors]

# Writing into file
file = open("geometry.txt", "w")
write(file, "x, y, z, size, front\n")
for i in 1:length(fx)
    fstring = string(fx[i]*d) * ", " * string(fy[i]*d) * ", " * string(zdif) * ", " * string(size) * ", true\n"
    write(file, fstring)
    bstring = string(bx[i]*d) * ", " * string(by[i]*d) * ", " * string(-zdif) * ", " * string(size) * ", false\n"
    write(file, bstring)
end
close(file)
