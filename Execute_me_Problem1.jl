# Pkg.add("LaTeXStrings")
using LaTeXStrings
using DifferentialEquations
using Plots
gr(show = true)

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout
AbstractPlotting.inline!(true)

# Model parameters
ν = 0.0001

#--------  Notch-Delta Two Cell Model -------------
#--------      Assuming ν << 1        -------------
function NotchDelta!(dd,d,ν,t)
    dd[1] = (1/(1+10*(d[2]^2/(0.1 + d[2]^2)^2)) - d[1])*ν
    dd[2] = (1/(1+10*(d[1]^2/(0.1 + d[1]^2)^2)) - d[2])*ν
end

# ------------- SOLVE THE MODEL WITH DIFFERENTIALEQUATIONS.jl -------------------
D₁ = 0.80                # initial value of Delta in cell 1
D₂ = 0.30                # initial value of Delta in cell 2
d0 = [D₁; D₂]       # initial state vector
tspan = (0.0,100000.0)                     #time interval (start time, end time)
prob = ODEProblem(NotchDelta!,d0,tspan, (ν))     #Create an ODE problem for the fxn
sol = solve(prob)                       #Solve the system

#--------------- PLOT REMPORAL RESULTS ------------------------------------------
#Plot the results; D₁ & D₂ vs time
plt1 = plot(sol, xaxis="time", yaxis = "D₁, D₂", label=["D₁" "D₂"])
display(plt1)

#--------  Notch-Delta Two Cell Dynamics -------------
#--------      Assuming ν << 1           -------------
function ND_dynamics(D₁,D₂)
    ν = 0.0001
    dD₁ = (1/(1+10*(D₂^2/(0.1 + D₂^2))^2) - D₁)*ν
    dD₂ = (1/(1+10*(D₁^2/(0.1 + D₁^2))^2) - D₂)*ν
    return Point(dD₁,dD₂)
end

# Construct the streamplot
f2 = Figure(backgroundcolor = RGBf0(0.98, 0.98, 0.98),resolution = (800, 800))
ax1 = f2[1,1] = Axis(f2)
ax1.ylabel = "D₁"
ax1.xlabel = "D₂"
str1 = streamplot!(ax1, ND_dynamics, 0..1, 0..1, colormap = :plasma, 
    gridsize= (50,50), arrow_size = 0.005)

# Add Nullclines
x1 = collect(LinRange(0,1,50))
y1 = 1 ./ (1 .+ 10*(x1.^2 ./ (0.1 .+ x1.^2)).^2)
line1 = lines!(ax1,x1,y1, color = (:green,0.9), linewidth = 1.5)

y2 = collect(LinRange(0,1,50))
x2 = 1 ./ (1 .+ 10*(y2.^2 ./ (0.1 .+ y2.^2)).^2)
line2 = lines!(ax1,x2,y2, color = (:blue,0.9), linewidth = 1.5)

# Display the plot
axislegend(ax1,[line1,line2],["dD₁/dτ=0","dD₂/dτ=0"],"Nullclines", position = :rc)
Makie.xlims!(ax1, [0.0, 1.0])
Makie.ylims!(ax1, [0.0, 1.0])
plt2 = f2.scene
ax1.title = "Two-Cell Delta Expression Phase Portrait"
display(f2)

# Save the plot
save("ND_Field.png", plt2)