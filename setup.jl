# Setup file for double pendulum project
# Run this once to install all required packages

using Pkg

# List of required packages
packages = [
    "OrdinaryDiffEq",  # Plus léger que DifferentialEquations
    "Plots",
    "Optim",
    "LightXML"
]

println("Installing required packages...")
for pkg in packages
    if !haskey(Pkg.project().dependencies, pkg)
        println("Installing $pkg...")
        Pkg.add(pkg)
    else
        println("$pkg already installed")
    end
end

println("\nAll packages installed successfully!")