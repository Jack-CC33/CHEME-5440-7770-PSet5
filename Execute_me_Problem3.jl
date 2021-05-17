using DifferentialEquations     # Include DifferentialEquations.jl
using Plots                     # Include Plots.jl for plotting
gr(show = true)  # Use the gr backend for plotting and show plots

# Model parameters
Vmax = 
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

function a1_plus(l)
    return 1/(1+l)
end

function a1_minus(l)
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
    du[3] = α1_plus(l)*u[2] - α1_minus(l)*u[3] + dbp*u[7] - adp*u[3]*u[5] + db*u[6] + ab*u[3]*u[4]
    du[4] = (db+kb+β1(l))*u[6] - (ab+k_plus)*u[3]*u[4] + k_minus*u[5]
    du[5] = (dbp+kbp+β1(l))*u[7] - abp*u[3]*u[5] - k_minus*u[5] + k_plus*u[3]*u[4]
    du[6] = ab*u[3]*u[4] - (db+kb+β1(l))*u[6]
    du[7] = abp*u[3]*u[5] - (dbp+kbp+β1(l))*u[7]
end

# ------------- SOLVE THE MODEL WITH DIFFERENTIALEQUATIONS.jl -------------------
U₁ = 1.0                # initial value of Eₒ
U₂ = 3.0                # initial value of E₁
U₃ = 0.0                # initial value of E₁*
U₄ = 0.0                # initial value of B
U₅ = 0.0                # initial value of Bₚ
u₆ = 0.0                # initial value of {E₁*B}
U₇ = 0.0                # initial value of {E₁*Bₚ}
u0 = [U₁; U₂; U₃; U₄; U₅; U₆; U₇]       # initial state vector
tspan = (0.0,100.0)                     #time interval (start time, end time)
prob = ODEProblem(Barkai!,u0,tspan)     #Create an ODE problem for the Lorenz fxn
sol = solve(prob)                       #Solve the system


#Plot the results; X, Y, and Z vs time
plt1 = plot(sol, xaxis="time", yaxis = "Concentrations", label=["Eₒ" "E₁" "E₁*" "B" "Bₚ" "{E₁*B}" "{E₁*Bₚ}"])
display(plt1)