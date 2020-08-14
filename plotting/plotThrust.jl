# Plot Thrust-Specific Constraints

using TrajOptSOCPs, Plots

function plotMagThrust_StartEnd(uListStart, uListEnd, maxThrust = 0)
    plt = plot()

    if maxThrust != 0
        hline!([maxThrust], linestyle = :dash, linecolor = :grey,
                            label = "Max Thrust")
    end

    colorStart = 1
    colorEnd = 2

    plot!(norm.(uListStart), label = "Initial Trajectory",
                        markershape = :square, markercolor = colorStart,
                        linecolor = colorStart)
    plot!(norm.(uListEnd), label = "Final Trajectory",
                        markershape = :diamond, markercolor = colorEnd,
                        linecolor = colorEnd)
    title!("Thrust Magnitude")
    ylabel!("Thrust")
    xlabel!("Time")

    # Set the lower limit to zero to be more accurate in representation
    ylims!(0.0, ylims()[2])
    xlims!(0.0, xlims()[2])

    return plt
end

function getAngle3D(u)
    if norm(u[3]) == 0.0
        # hit the tan singularity
        return pi / 2
    else
        return atan(norm([u[1]; u[2]]), -u[3])
    end
end

function plotAngle3D_StartEnd(uListStart, uListEnd, maxAngle = 0, deg = true)
    plt = plot()

    if maxAngle != 0
        hline!([maxAngle, -maxAngle], linestyle = :dash, linecolor = :grey,
                            label = "Max Angle")
    end

    colorStart = 1
    colorEnd = 2

    if deg
        # get the rad to deg conversion
        conversion = 180 / pi
    else
        # leave as rad
        conversion = 1
    end


    plot!(getAngle3D.(uListStart) * conversion, label = "Initial Trajectory",
                        markershape = :square, markercolor = colorStart,
                        linecolor = colorStart)
    plot!(getAngle3D.(uListEnd) * conversion, label = "Final Trajectory",
                        markershape = :diamond, markercolor = colorEnd,
                        linecolor = colorEnd)
    title!("Thrust Vector Angle")
    if deg
        ylabel!("Angle (deg)")
    else
        ylabel!("Angle (rad)")
    end
    xlabel!("Time")

    # Set the lower limit to zero to be more accurate in representation
    xlims!(0.0, xlims()[2])

    return plt
end
