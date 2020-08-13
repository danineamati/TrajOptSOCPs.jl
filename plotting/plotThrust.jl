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
