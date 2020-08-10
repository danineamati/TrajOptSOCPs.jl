# Trajectory Parsing

"""
    parseTrajectory(sList, vList, uList)

Stores a parsed trajectory by separating the positions (sList), the velocities
(vList), and the controls (uList).

See [`getParseTrajectory`](@ref) for more info.
"""
struct parseTrajectory
    sList
    vList
    uList
end

"""
    getParseTrajectory(traj, nDim)

Converts a trajectory of the form `[x0; u0; x1; u1; ...]` where `xk = [sk; vk]`
to a [`parseTrajectory`](@ref) struct.
"""
function getParseTrajectory(traj, nDim)
    sList = []
    vList = []
    uList = []

    kStart = 1
    maxK = length(traj)
    # println("Max = $maxK")

    while kStart < (maxK - 2 * nDim)
        # println("Accessing $kStart to $(kStart + 3 * nDim - 1)")
        push!(sList, traj[kStart:kStart + nDim - 1])
        push!(vList, traj[kStart + nDim:kStart + 2 * nDim - 1])
        push!(uList, traj[kStart + 2 * nDim:kStart + 3 * nDim - 1])

        kStart = kStart + 3 * nDim
    end

    # Add the final state (has no control step)
    # println("Last Access is $(kStart + 2 * nDim - 1)")
    push!(sList, traj[kStart:kStart + nDim - 1])
    push!(vList, traj[kStart + nDim:kStart + 2 * nDim - 1])

    return parseTrajectory(sList, vList, uList)
end

"""
    splitDimensions(stateList)

Converts a list of `n` dimensional states into `n` lists of 1 dimension.
"""
function splitDimensions(stateList)
    nDim = size(stateList[1], 1) # How many dimensions?

    splitStateList = []

    for d in 1:nDim
        push!(splitStateList, [x[d] for x in stateList])
    end

    return splitStateList
end
