# Generates a batch of plots

using TrajOptSOCPs

include("plotConstraintViolation.jl")
include("plotObjective.jl")
include("plotTrajectory.jl")
include("plotTrajectory3D.jl")
include("plotThrust.jl")


function batchPlot(trajStates, cM::constraintManager, nDim::Int64,
                            penalty::Float64 = 1.0, maxThrust::Float64 = 0.0,
                            maxAngleDeg::Float64 = 0.0)
    ptList = [getParseTrajectory(traj, nDim) for traj in trajStates]
    # println("ptList: ")
    # display(ptList)

    ptrajStart = ptList[1]
    ptrajLast = ptList[end]

    if nDim == 2
        pltTraj = plotTrajPos2D_Multiple(ptList)
        display(pltTraj)
    elseif nDim == 3
        pltTraj = plotTrajPos3D_Multiple(ptList)
        display(pltTraj)
    end


    pltCV, pltCV2 = plotConstraintViolation(cM, trajStates, penalty)
    display(pltCV)
    display(pltCV2)

    pltObj = plotObjective(costFun, trajStates)
    display(pltObj)

    # plotSVUTime_Simple(ptList[1])
    if nDim == 2
        plts, pltv, pltu = plotSVUTime_StartEnd([ptrajStart, ptrajLast])
    elseif nDim == 3
        plts, pltv, pltu = plotSVUTime3D_StartEnd([ptrajStart, ptrajLast])
    end

    pltMagT = plotMagThrust_StartEnd(ptrajStart.uList, ptrajLast.uList,
                                     maxThrust)
    display(pltMagT)

    pltAngleT = plotAngle3D_StartEnd(ptrajStart.uList, ptrajLast.uList,
                                     maxAngleDeg)
    display(pltAngleT)


    return pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu, pltMagT, pltAngleT
end


function saveBulk(pltTraj, pltCV, pltCV2, pltObj, plts, pltv, pltu,
                    pltMagT, pltAngleT, header::String)
    savefig(pltTraj, header * "Trajectory")
    savefig(pltCV, header * "ConstraintViolation")
    savefig(pltCV2, header * "DynamicsViolation")
    savefig(pltObj, header * "Objective")
    savefig(plts, header * "PositionTime")
    savefig(pltv, header * "VelocityTime")
    savefig(pltu, header * "ControlsTime")
    savefig(pltMagT, header * "ThrustMagTime")
    savefig(pltAngleT, header * "ThrustAngleTime")
end
