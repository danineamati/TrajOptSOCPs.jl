# Similar to plotTrajectory.jl, but in 3D

using Plots

function plotTrajPos3D_Multiple(ptList, showN = 10)
    plt = plot3d(legend = :outerright, zlabel = "Z (km)")

    if length(ptList) >= showN
        indices = floor.(Int, LinRange(1, length(ptList), showN))
    else
        indices = collect(1:length(ptList))
    end

    println("Plotting the following indices: ")
    println(indices)

    minXY = 0
    maxXY = 0

    for ind in indices
        pt = ptList[ind]
        xyzList = splitDimensions(pt.sList)
        xList = xyzList[1]
        yList = xyzList[2]
        zList = xyzList[3]
        plot3d!(xList, yList, zList, markershape = :circle,
                                label = "Traj $ind")

        minXY = min(minXY, minimum(xList), minimum(yList))
        maxXY = max(maxXY, maximum(xList), maximum(yList))
    end

    xlabel!("X (km)")
    ylabel!("Y (km)")
    # zlabel!("Z (km)")
    title!("Trajectory Convergence")
    # yflip!()

    # Check that the XY plane looks square
    xlims!(minXY, maxXY)
    ylims!(minXY, maxXY)

    return plt
end

function plotSVUTime3D_StartEnd(ptlist)
    return plot(), plot(), plot()
end
