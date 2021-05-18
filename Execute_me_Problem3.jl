using DifferentialEquations     # Include DifferentialEquations.jl
using Plots                     # Include Plots.jl for plotting
using Printf
gr(show = true)                 # Use the gr backend for plotting and show plots

# Model parameters (all units in terms of uM^n and s^m)
Vmax = 0.13  # Estimated value
ab = 1000
db = 1
kb = 0
abp = 100
dbp = 0.01
kbp = 1
k_plus = 1
k_minus = 1
α0_plus = 10
α0_minus = 0

function β1(l)
    return 2.5*l/(1+l)
end

function α1_plus(l)
    return 1/(1+l)
end

function α1_minus(l)
    return l/(1+l)
end

E_conc = 10
R_conc = 0.2
B_conc = 2

# -------------------------- Barkai Model -------------------------------------
# du: Diffrerential equations, [dEo/dt, dE1/dt, dE1*/dt, dB/dt, dBp/dt, d{E1*B}/dt, d{E1*Bp}/dt]
# u: Time-dependent variables, [Eo, E1, E1*, B, Bp, {E1*B}, {E1*Bp}]
# p: Additional model parameters (none)
# t: time 
function Barkai!(du,u,l,t)
    du[1] = kbp*u[7] + kb*u[6] - Vmax
    du[2] = Vmax + α1_minus(l)*u[3] - α1_plus(l)*u[2] + β1(l)*(u[6]+u[7])
    du[3] = α1_plus(l)*u[2] - α1_minus(l)*u[3] + dbp*u[7] - abp*u[3]*u[5] + db*u[6] - ab*u[3]*u[4]
    du[4] = (db+kb+β1(l))*u[6] - (ab+k_plus)*u[3]*u[4] + k_minus*u[5]
    du[5] = (dbp+kbp+β1(l))*u[7] - abp*u[3]*u[5] - k_minus*u[5] + k_plus*u[3]*u[4]
    du[6] = ab*u[3]*u[4] - (db+kb+β1(l))*u[6]
    du[7] = abp*u[3]*u[5] - (dbp+kbp+β1(l))*u[7]
end

# ------------- SOLVE THE MODEL WITH DIFFERENTIALEQUATIONS.jl -------------------
U₁ = 7.170677164815677                      # initial value of Eₒ
U₂ = 0.13                                   # initial value of E₁
U₃ = 0.7038387614777921                     # initial value of E₁*
U₄ = 0.0026504423680775635                  # initial value of B
U₅ = 0.0018654840737143306                  # initial value of Bₚ
U₆ = 1.8654840737159795                     # initial value of {E₁*B}
U₇ = 0.1299999999999626                     # initial value of {E₁*Bₚ}
u0 = [U₁; U₂; U₃; U₄; U₅; U₆; U₇]           # initial state vector
tspan = (0.0,100.0)                       #time interval (start time, end time)
prob = ODEProblem(Barkai!,u0,tspan,0.01)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                           #Solve the system


#Plot the results; X, Y, and Z vs time
h = []
for i ∈ 1:size(sol.u)[1]
    append!(h, [sol.u[i][3]/sol.u[1][3]])
end
plt1 = plot(sol, xaxis="time", yaxis = "Concentrations", label=["Eₒ" "E₁" "E₁*" "B" "Bₚ" "{E₁*B}" "{E₁*Bₚ}"])
display(plt1)

# L = 0.0 → 0.1 μM
plt2 = plot(sol.t,h, xaxis="time", yaxis="Activity")
title!(plt2,"[L] = 0.0 → 0.01 μM")
ylims!(plt2,(0.0,1.3))
display(plt2)


# L = 0.0 → 0.1 μM
tspan = (0.0,300.0)                       #time interval (start time, end time)
prob = ODEProblem(Barkai!,u0,tspan,0.1)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                           #Solve the system
h = []
for i ∈ 1:size(sol.u)[1]
    append!(h, [sol.u[i][3]/sol.u[1][3]])
end
plt3 = plot(sol.t,h, xaxis="time", yaxis="Activity")
title!(plt3,"[L] = 0.0 → 0.1 μM")
ylims!(plt3,(0.0,2.0))
display(plt3)


# L = 0 → 1 μM
tspan = (0.0,400.0)                       #time interval (start time, end time)
prob = ODEProblem(Barkai!,u0,tspan,1)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                           #Solve the system
h = []
for i ∈ 1:size(sol.u)[1]
    append!(h, [sol.u[i][3]/sol.u[1][3]])
end
plt4 = plot(sol.t,h, xaxis="time", yaxis="Activity")
title!(plt4,"[L] = 0 → 1 μM")
ylims!(plt4,(0.0,2.0))
display(plt4)


# L = 0 → 10 μM
tspan = (0.0,1000.0)                       #time interval (start time, end time)
prob = ODEProblem(Barkai!,u0,tspan,10)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                           #Solve the system
h = []
for i ∈ 1:size(sol.u)[1]
    append!(h, [sol.u[i][3]/sol.u[1][3]])
end
plt5 = plot(sol.t,h, xaxis="time", yaxis="Activity")
title!(plt5,"[L] = 0 → 10 μM")
ylims!(plt5,(0.0,2.0))
display(plt5)



# Print final state
#U_final =  sol.u[end]
#print(U_final)

# save figures
savefig(plt1, "./Problem3_10nM_Ligand_Full_Response.png") 
savefig(plt2, "./Problem3_10nM_Ligand_Activity.png") 
savefig(plt3, "./Problem3_100nM_Ligand_Activity.png") 
savefig(plt4, "./Problem3_1uM_Ligand_Activity.png") 
savefig(plt5, "./Problem3_10uM_Ligand_Activity.png") 