# Creates a ground constraint

using SparseArrays

include("..\\constraints\\constraints.jl")


"""
    makeGroundConstraint(NSteps::Int64, nDim::Int64, index = 2)

Creates an [`AL_AffineInequality`](@ref) constraint such that the rocket stays
above the ground. (Constraint has form Ax - b ≤ 0). Ignores the last point
in the trajectory.

Usage: `NSteps` states the number of timesteps in the trajectory. `nDim` states
the dimensionality of the problem (e.g. 1D, 2D, or 3D). Index describes which
direction is up (e.g. index = 2 means "y" is up).
"""
function makeGroundConstraint(NSteps::Int64, nDim::Int64, index = 2)
    # Safegaurd to assue the input is reasonable
    if index > nDim
        print("Invalid index. Changed from $index to ")
        index = min(nDim, index)
        println(index)
    end

    # Build up the sparse matrices
    sizeTraj = 3 * NSteps * nDim + 2 * nDim
    AMat = spzeros(NSteps, sizeTraj)
    bVec = spzeros(NSteps)

    # AMat = zeros(NSteps, sizeTraj)
    # bVec = zeros(NSteps)

    # Iterate through and select each index that needs to be flipped to -1
    currInd = index
    currRow = 1
    while currRow <= NSteps
        AMat[currRow, currInd] = -1

        currRow += 1
        currInd += 3 * nDim
    end

    # return  the result as a constraint struct
    return AL_AffineInequality(AMat, bVec)

end
