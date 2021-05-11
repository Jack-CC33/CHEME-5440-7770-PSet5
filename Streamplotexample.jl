# Example for constructing a streamplot in Julia for phase portraits
# CHEME5440/7770 - Spring 2021
# The example makes use of Makie.jl and AbstractPlotting.jl for plotting
# To install the packages, use the following commands inside the Julia REPL:
# using Pkg
# Pkg.add("Makie")
# Pkg.add("AbstractPlotting")

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout
AbstractPlotting.inline!(true)

# Model for precise adaptation
# Xa: active receptor state
# Xi: inactive receptor state
function precise_adapt(Xa, Xi)
    Vrmax = 1               #Vmax for conversion of X to Xa by R (kcat,r*R)
    Vbmax = 4               #Vmax for conversion of Xa to X (kcat,b*B)
    Ka = 1                  #Michaelis-menten constant for CheB
    kon = 1                 #rate constant for ligand binding
    koff = 1                #rate constant for ligand unbinding
    c = 4                   #ligand concentration
    
    u = Vrmax - Vbmax*Xa/(Ka + Xa) - kon*c*Xa + koff*Xi     #dXa/dt
    v = kon*c*Xa - koff*Xi                                  #dXi/dt
    
    return Point(u,v)
end

# Construct the streamplot
plt1 = Scene(resolution =(400,400))
streamplot!(plt1, precise_adapt, 0..4, 0..4, colormap = :plasma, 
    gridsize= (32,32), arrow_size = 0.05)

# Display the plot
display(plt1)

# Save the plot
save("odeField.png", plt1)