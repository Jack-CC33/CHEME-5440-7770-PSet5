# Example for use of DifferentialEquations.jl
# CHEME5440/7770 - Spring 2021
# To install the package, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("DifferentialEquations")
# For this example, you will also need Plots; To add:
# Pkg.add("Plots")
# ------------------------------------------------------------ #
# Lorenz system example
# dX/dt = σ(Y-X)
# dY/dt = X(ρ - Z) - Y
# dZ/dt = XY - βZ
# Based on DifferentialEquations.jl tutorial

using DifferentialEquations     # Include DifferentialEquations.jl
using Plots                     # Include Plots.jl for plotting
gr(show = true)  # Use the gr backend for plotting and show plots

# Model parameters
σ = 5.5
ρ = 22.0   #.0 indicates to the compiler that ρ is a floating-point number (i.e. not an integer)
β = 7/3

# -------------------------- Lorenz Model -------------------------------------
# du: Diffrerential equations, [dX/dt, dY/dt, dZ/dt]
# u: Time-dependent variables, [X, Y, Z]; u[1] = xX; u[2] = Y; u[3] = Z
# p: Additional model parameters (none in this example)
# t: time 
# Note "!" point after function name is a Julia convention that indicates 
# that the function will modify values in one or more of the function arguments.  In this case,
# lorenz! will modify values in the input vector, du
function lorenz!(du,u,p,t)
 du[1] = σ*(u[2]-u[1])                  #dX/dt = σ(Y-X)
 du[2] = u[1]*(ρ-u[3]) - u[2]           #dY/dt = X(ρ - Z) - Y
 du[3] = u[1]*u[2] - β*u[3]             #dZ/dt = XY - βZ
end

# ------------- SOLVE THE MODEL WITH DIFFERENTIALEQUATIONS.jl -------------------
X₀ = 1.0                # initial value of X
Y₀ = 3.0                # initial value of Y
Z₀ = 0.0                # initial value of Z
u0 = [X₀; Y₀; Z₀]       # initial state vector
tspan = (0.0,100.0)                     #time interval (start time, end time)
prob = ODEProblem(lorenz!,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                       #Solve the system

# ------------- MAKE SOME PLOTS OF THE RESULTS WITH PLOTS.jl ---------------------
#Plot the results; X, Y, and Z vs time
plt1 = plot(sol, xaxis="time", yaxis = "X, Y, or Z", label=["X" "Y" "Z"])
display(plt1)
#Plot the results; the vars=(0,3) argument specifies to plot Z (column 3 of sol)
#vs t
plt2 = plot(sol,vars=(0,3), xaxis="time", yaxis = "Z", legend = false)
display(plt2)
#Plot the results; the vars=(1,2,3) argument specifies to plot X vs Y vs Z
plt3 = plot(sol,vars=(1,2,3), xaxis="X", yaxis="Y", zaxis="Z", legend = false)
display(plt3)
#Save the three plots as PNG files
savefig(plt1, "./plot1.png")    
savefig(plt2, "./plot2.png")  
savefig(plt3, "./plot3.png")    