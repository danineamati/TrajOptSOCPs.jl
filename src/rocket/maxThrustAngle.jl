# Add Max Trust Angle Constraints

using SparseArrays, LinearAlgebra

include("../constraints/constraints.jl")


function getAlpha(theta, deg = true)
    if deg
        return -tand(theta)
    end
    return -tan(theta)
end


"""
    makeMaxThrustConstraint(NSteps::Int64, nDim::Int64, thetaMax::Float64)

Creates an [`AL_Multiple_simpleCone`](@ref) constraint such that the trust
(or equivalently control) vector angle is capped by `thetaMax`.

Usage: `NSteps` states the number of timesteps in the trajectory. `nDim` states
the dimensionality of the problem (e.g. 1D, 2D, or 3D). `thetaMax` is the max
thrust angle.
"""
function makeMaxAngleConstraint(NSteps::Int64, nDim::Int64, thetaMax::Float64,
                                    deg = true)
    # sizeTraj = 3 * NSteps * nDim + 2 * nDim

    # This helper function gets the index of ux for any k
    startUInd(k) = (2 * nDim + 1) + (3 * nDim * (k - 1))

    indicatorList = Int64[startUInd(k) for k in 1:NSteps]

    sInd = Diagonal([ones(Int64, nDim - 1); 0]) # e.g. Diagonal([1; 1; 0])
    tInd = [zeros(Int64, nDim - 1); 1]          # e.g. [0; 0; 1]

    alpha = getAlpha(thetaMax, deg)

    return makeAL_Multiple_AngleCone(alpha, sInd, tInd,
                                        nDim, indicatorList)
end


function makeMaxAngleConstraint(NSteps::Int64, nDim::Int64, thetaMax::Float64,
                                    sIndicate::Diagonal{Int64,Array{Int64,1}},
                                    tIndicate::Array{Int64,1}, deg = true)
    # sizeTraj = 3 * NSteps * nDim + 2 * nDim

    # This helper function gets the index of ux for any k
    startUInd(k) = (2 * nDim + 1) + (3 * nDim * (k - 1))

    indicatorList = Int64[startUInd(k) for k in 1:NSteps]

    alpha = getAlpha(thetaMax, deg)

    return makeAL_Multiple_simpleCone(alpha, sIndicate, tIndicate,
                                        nDim, indicatorList)
end
