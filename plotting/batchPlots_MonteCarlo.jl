# Generates a batch of plots for a MonteCarlo run.

using TrajOptSOCPs

include("plotTrajectory.jl")
include("plotTrajectory3D.jl")
include("plotConstraintViolation_MonteCarlo.jl")


function batchPlotMonteCarlo(trajList, dynVio, constVio, nDim::Int64,
                            penalty::Float64 = 1.0, maxThrust::Float64 = 0.0,
                            maxAngleDeg::Float64 = 0.0)

    # First and foremost, we want to plot the trajectories
    ptList = [getParseTrajectory(traj, nDim) for traj in trajList]

    if nDim == 2
        pltTraj = plotTrajPos2D_Multiple(ptList)
        display(pltTraj)
    elseif nDim == 3
        pltTraj = plotTrajPos3D_Multiple(ptList)
        display(pltTraj)
    end

    pltVio = plotDynConstVio(dynVio, constVio, true)

    pltList = [pltTraj, pltVio]
    pltLabels = ["Trajectory", "Violation"]

    # Return the result as a dictionary
    return Dict(zip(pltLabels, pltList))
end
