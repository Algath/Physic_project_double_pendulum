using HypertextLiterals
using PlutoUI

md"""
# Double pendulum simulation

This is a simulation of a double pendulum using Julia.

The objectif is too recreate a double pendulum simulation based on videos and to predict his mouvement during 2 seconds past the video end.
"""

md"""
## Theory
A double pendulum is also know as Chaotic pendulum because of its sensibility to initial conditions.

The particularity of this system is his behavior. It can vary from simple periodic motion to a more complex and chaotic one with a small change of initial angles or velocities.

"""

md"""
## length of the pendulums (via gimps)
pendulum one = 91.74 mm
pendulum two = 69.33 mm
"""

md"""
## Parameters
- length of the first pendulum: L1 = 91.74 mm
- length of the second pendulum: L2 =  69.33 mm
- angles at t=0: θ1 = 181.5° = 3.159, θ2 = 183.1° = 3.194 (θ = 0° is the vertical down position)
"""

function degrees_to_radians(degrees)
    return degrees * (π / 180)
end
function Lagrangien(θ1_0, θ2_0, m1, m2, g, r1, r2)
end
function position(r1, θ1_0, L1, r2, θ2_0, L2)
    r1[1] = L1*sin(θ1_0)
    r1[2] = -L1*cos(θ1_0)

    r2[1] = r1[1] + L2*sin(θ2_0)
    r2[2] = r1[2] - L2*cos(θ2_0)
end

L1 = 91.74
L2 = 69.33
θ1_0 = degrees_to_radians(181.5)
θ2_0 = degrees_to_radians(183.1)
g = 9.81  # Acceleration due to gravity (m/s^2)
m1 = 0
m2 = 0

r1 = (0, 0)
r2 = (0, 0)