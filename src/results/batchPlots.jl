# Generates a batch of plots

include("plotTrajectory.jl")
include("plotConstraintViolation.jl")
include("plotObjective.jl")
include("trajectoryParsing.jl")


function batchPlot(trajStates, cM::constraintManager, nDim::Int64,
                            penalty::Float64 = 1.0)
    ptList = [getParseTrajectory(traj, nDim) for traj in trajStates]
    # println("ptList: ")
    # display(ptList)

    ptrajStart = ptList[1]
    ptrajLast = ptList[end]

    pltTraj = plotTrajPos2D_Multiple(ptList)
    display(pltTraj)

    pltCV, pltCV2 = plotConstraintViolation(cM, trajStates, penalty)
    display(pltCV)
    display(pltCV2)

    pltObj = plotObjective(costFun, trajStates)
    display(pltObj)

    # plotSVUTime_Simple(ptList[1])
    plts, pltv, pltu = plotSVUTime_StartEnd([ptrajStart, ptrajLast])

    return pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu
end


function saveBulk(pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu,
                    header::String)
    savefig(pltTraj, header * "Trajectory")
    savefig(pltCV, header * "ConstraintViolation")
    savefig(pltCV2, header * "DynamicsViolation")
    savefig(pltObj, header * "Objective")
    savefig(plts, header * "PositionTime")
    savefig(pltv, header * "VelocityTime")
    savefig(pltu, header * "ControlsTime")
end
