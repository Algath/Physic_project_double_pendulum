# Import all required packages for the double pendulum simulation
# Include this file in your main script with: include("imports.jl")
include("setup.jl")

using OrdinaryDiffEq         # Pour résoudre les EDO (plus léger que DifferentialEquations)
using Plots                  # Pour la visualisation et animations
using Optim                  # Pour l'optimisation des masses
using Statistics             # Pour mean() dans RMSE/R² (standard library)
using LightXML               # Pour parser les fichiers Tracker .trk (XML)

println("All packages loaded successfully!")
