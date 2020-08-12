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

function plotSVUTime3D_Multiple(sListStack, vListStack, uListStack)

    # Variables needed for later
    nDim = size(sListStack[1][1], 1)
    sxyzLabels = ["x", "y", "z"]
    vxyzLabels = ["vx", "vy", "vz"]
    uxyzLabels = ["ux", "uy", "uz"]

    seLabels = ["Initial", "Final"]

    markerShapeArr = [:rect, :diamond]
    markerSizeArr = [8; 4]
    markerColorArr = [:blue, :purple, :darkgreen, :orange2, :gold, :pink]

    # Plot the Positions
    plt_s = plot()

    colorSelect = 1
    for (ind, sList) in enumerate(sListStack)
        xyzList = splitDimensions(sList)
        for d in 1:nDim
            thisLabel = sxyzLabels[d] * " for " * seLabels[ind] * " Traj"
            plot!(xyzList[d], markershape = markerShapeArr[ind],
                     markersize = markerSizeArr[ind],
                     markercolor = markerColorArr[colorSelect],
                     linecolor = markerColorArr[colorSelect],
                     label = thisLabel)
            colorSelect += 1
        end
    end
    title!("Change in position over time")
    xlabel!("Time")
    ylabel!("Position")
    display(plt_s)

    # Plot the Velocities
    plt_v = plot()
    colorSelect = 1
    for (ind, vList) in enumerate(vListStack)
        vxyzList = splitDimensions(vList)
        for d in 1:nDim
            thisLabel = vxyzLabels[d] * " for " * seLabels[ind] * " Traj"
            plot!(vxyzList[d], markershape = markerShapeArr[ind],
                     markersize = markerSizeArr[ind],
                     markercolor = markerColorArr[colorSelect],
                     linecolor = markerColorArr[colorSelect],
                     label = thisLabel)
            colorSelect += 1
        end
    end
    title!("Change in velocity over time")
    xlabel!("Time")
    ylabel!("Velocity")
    display(plt_v)

    # Plot the Controls
    colorSelect = 1
    markerColorArrU = [ markerColorArr[1], markerColorArr[4],
                        markerColorArr[2], markerColorArr[5],
                        markerColorArr[3], markerColorArr[6]]

    uxyzListSplit = splitDimensions.(uListStack)

    uPltList = []

    for d in 1:nDim
        pltNew = plot(legend = :outerright)
        for ind in 1:length(uxyzListSplit)
            thisLabel = uxyzLabels[d] * " for " * seLabels[ind] * " Traj"
            uxyzList = uxyzListSplit[ind][d]
            plot!(uxyzList, markershape = markerShapeArr[ind],
                     markersize = markerSizeArr[ind],
                     markercolor = markerColorArrU[colorSelect],
                     linecolor = markerColorArrU[colorSelect],
                     label = thisLabel)
            colorSelect += 1
        end
        push!(uPltList, pltNew)
        # display(pltNew)
    end

    plt_u = plot(uPltList..., layout = (length(uPltList), 1))

    title!("Change in thrust over time")
    xlabel!("Time")
    ylabel!("Contol (Thrust)")
    display(plt_u)

    return  plt_s, plt_v, plt_u
end

function plotSVUTime3D_StartEnd(ptList::Array{parseTrajectory, 1})
    ptStart = ptList[1]
    ptEnd = ptList[end]

    sListStack = [ptStart.sList, ptEnd.sList]
    vListStack = [ptStart.vList, ptEnd.vList]
    uListStack = [ptStart.uList, ptEnd.uList]

    plotSVUTime3D_Multiple(sListStack, vListStack, uListStack)
end
