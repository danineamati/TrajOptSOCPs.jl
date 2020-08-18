# Animate a 3D Trajectory

getbounds(list) = [minimum(list), maximum(list)]

function animTraj3D(pt::parseTrajectory, lag = 5, land = [0; 0; 0])
    xyzList = splitDimensions(pt.sList)

    x = xyzList[1]
    y = xyzList[2]
    z = xyzList[3]

    x_b = getbounds(x)
    y_b = getbounds(y)
    z_b = getbounds(z)

    cube = getbounds([x_b..., y_b...])

    N = length(x)

    my_anim = @animate for i in 1:(N + lag)
        j = min(i, N)
        plot3d(x[1:j], y[1:j], z[1:j], markershape = :circle,
                        zlabel = "Z (m)", ylabel = "Y (m)", xlabel = "X(m)",
                        legend = :none, zlims = (z_b[1], z_b[2]))
        scatter3d!([land[1]], [land[2]], [land[3]],
                        markershape = :xcross, markercolor = :red,
                        markersize = 8)
        title!("Trajectory over time")

        xlims!(cube...)
        ylims!(cube...)
    end

    return my_anim
end
