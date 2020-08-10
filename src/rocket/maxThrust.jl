# Add Max Trust Constraints

using SparseArrays

include("..\\constraints\\constraints.jl")

"""
    makeMaxThrustConstraint(NSteps::Int64, nDim::Int64, tMax::Float64)

Creates an [`AL_Multiple_simpleCone`](@ref) constraint such that the trust
(or equivalently control) vector is capped by `tMax`.

Usage: `NSteps` states the number of timesteps in the trajectory. `nDim` states
the dimensionality of the problem (e.g. 1D, 2D, or 3D). `tMax` is the max
thrust.
"""
function makeMaxThrustConstraint(NSteps::Int64, nDim::Int64, tMax::Float64)
    # sizeTraj = 3 * NSteps * nDim + 2 * nDim

    # This helper function gets the index of ux for any k
    startUInd(k) = (2 * nDim + 1) + (3 * nDim * (k - 1))

    indicatorList = Int64[startUInd(k) for k in 1:NSteps]

    return makeAL_Multiple_simpleCone(tMax, nDim, indicatorList)
end
