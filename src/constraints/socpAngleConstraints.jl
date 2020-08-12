# SOCP Constraints for angle constraints

#=

These are constraints of the form
||[x1; x2]|| ≤ α x3
Since this is a trajectory, we also likely want to use this constraint
multiple times (e.g. at each step in the discretization)

This file handles these constraints such that it can be integrated into
the "constraintManager" framework.


Mathematics:
For the rocket soft-landing problem, we want to limit the angle that
engine can swing. We make the following assumptions:
1. The angle is usually to a small angle < 90°.
This holds because we don't want the engine gimbal to move too far from its
nominal axial position (i.e. engine exhause aligned with the rocket's body).
We also don't want the rocket control to go unstable due to a large torque
from a large angle.
2. The thrust vector is assumed to be _exactly_ in line with the engine
pointing vector.
This is a good approximation for a good rocket. There may be small variations,
but the engine is designed to direct the exhaust.
3. The rocket is oriented with the ground. (i.e. the long axial direction of
the rocket is perpendicular to the ground)
This is a good approximation near landing since the rocket landing gear must
downward facing to interface with the ground.

Then, we calculate the radial distance as ||ux|| in 2D and ||[ux; uy]|| in 3D.
The angle against the z component is
θ = arctan(||uxy|| / -uz)
The "-uz" is to have "-uz > 0" since the engine points downward. We want
θ < θ_max
Equivalently,
arctan(||uxy|| / -uz) < θ_max
||uxy|| < - tan(θ_max) uz

Call tan(θ_max) = α. Thus, we recover
||ux|| < - α uy in 2D
||[ux; uy]|| < - α uz in 3D

as desired.

=#


using SparseArrays

include("projections.jl")
