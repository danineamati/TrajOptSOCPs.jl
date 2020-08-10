# Plot Trajectories

include("trajectoryParsing.jl")


using Plots
pyplot()

function plotTrajPos2D_Simple(xList, yList)
    plt = plot(xList, yList, markershape = :circle)
    return plt
end

function plotTrajPos2D_Simple(pt::parseTrajectory)
    xyList = splitDimensions(pt.sList)
    xList = xyList[1]
    yList = xyList[2]
    return plotTrajPos2D_Simple(xList, yList)
end

function plotTrajPos2D_Multiple(ptList, showN = 10)
    plt = plot(legend = :outerright)

    if length(ptList) >= showN
        indices = floor.(Int, LinRange(1, length(ptList), showN))
    else
        indices = collect(1:length(ptList))
    end

    println("Plotting the following indices: ")
    println(indices)

    for ind in indices
        pt = ptList[ind]
        xyList = splitDimensions(pt.sList)
        xList = xyList[1]
        yList = xyList[2]
        plot!(xList, yList, markershape = :circle, label = "Traj $ind")
    end

    xlabel!("X (km)")
    ylabel!("Y (km)")
    title!("Test Trajectory")
    # yflip!()

    return plt
end

function plotSVUTime_Simple(sList, vList, uList)

    # Variables needed for later
    nDim = size(sList[1], 1)
    sxyzLabels = ["x", "y", "z"]
    vxyzLabels = ["vx", "vy", "vz"]
    uxyzLabels = ["ux", "uy", "uz"]

    # Plot the Positions
    plt_s = plot()
    xyzList = splitDimensions(sList)

    for d in 1:nDim
        plot!(xyzList[d], markershape = :circle, label = sxyzLabels[d])
    end
    title!("Change in position over time")
    xlabel!("Time")
    ylabel!("Position")
    display(plt_s)

    # Plot the Velocities
    plt_v = plot()
    vxyzList = splitDimensions(vList)

    for d in 1:nDim
        plot!(vxyzList[d], markershape = :circle, label = vxyzLabels[d])
    end
    title!("Change in velocity over time")
    xlabel!("Time")
    ylabel!("Velocity")
    display(plt_v)

    # Plot the Controls
    plt_u = plot()
    uxyzList = splitDimensions(uList)

    for d in 1:nDim
        plot!(uxyzList[d], markershape = :circle, label = uxyzLabels[d])
    end
    title!("Change in thrust over time")
    xlabel!("Time")
    ylabel!("Contol (Thrust)")
    display(plt_u)

    return  plt_s, plt_v, plt_u
end

function plotSVUTime_Simple(pt::parseTrajectory)
    return plotSVUTime_Simple(pt.sList, pt.vList, pt.uList)
end

function plotSVUTime_Multiple(sListStack, vListStack, uListStack)

    # Variables needed for later
    nDim = size(sListStack[1][1], 1)
    sxyzLabels = ["x", "y", "z"]
    vxyzLabels = ["vx", "vy", "vz"]
    uxyzLabels = ["ux", "uy", "uz"]

    markerShapeArr = [:rect, :diamond]
    markerSizeArr = [8; 4]
    markerColorArr = [:blue, :purple, :orange2, :gold]

    # Plot the Positions
    plt_s = plot()

    colorSelect = 1
    for (ind, sList) in enumerate(sListStack)
        xyzList = splitDimensions(sList)
        for d in 1:nDim
            thisLabel = sxyzLabels[d] * " for traj $ind"
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
            plot!(vxyzList[d], markershape = markerShapeArr[ind],
                     markersize = markerSizeArr[ind],
                     markercolor = markerColorArr[colorSelect],
                     linecolor = markerColorArr[colorSelect],
                     label = vxyzLabels[d])
            colorSelect += 1
        end
    end
    title!("Change in velocity over time")
    xlabel!("Time")
    ylabel!("Velocity")
    display(plt_v)

    # Plot the Controls
    colorSelect = 1
    markerColorArr = [:blue, :orange2, :purple, :gold]

    uxyzListSplit = splitDimensions.(uListStack)

    uPltList = []

    for d in 1:nDim
        pltNew = plot()
        for ind in 1:length(uxyzListSplit)
            uxyzList = uxyzListSplit[ind][d]
            plot!(uxyzList, markershape = markerShapeArr[ind],
                     markersize = markerSizeArr[ind],
                     markercolor = markerColorArr[colorSelect],
                     linecolor = markerColorArr[colorSelect],
                     label = uxyzLabels[d])
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

function plotSVUTime_StartEnd(ptList::Array{parseTrajectory, 1})
    ptStart = ptList[1]
    ptEnd = ptList[end]

    sListStack = [ptStart.sList, ptEnd.sList]
    vListStack = [ptStart.vList, ptEnd.vList]
    uListStack = [ptStart.uList, ptEnd.uList]

    plotSVUTime_Multiple(sListStack, vListStack, uListStack)
end
