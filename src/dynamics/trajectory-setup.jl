# Trajectory Set-Up




"""
    initializeTraj(x0::Array{Float64, 1}, xN::Array{Float64, 1}, NSteps::Int64)

Initializes the full `XU` vector (see example below) using a linear
interpolation from the initial point `x0` to the final point `xN`.

Take `x0 = [s0; v0]` and `xN = [sN; vN`] to be of size (`2n×1`). The full `XU`
vector size is `(3Nn + 2n)×1`.

For example, take a 2D system: `x0 = [sx0; sy0; vx0; vy0]` and
                               `xN = [sxN; syN; vxN, vyN]`.

The controls are `uk = [uxk; uyk]`. Thus, for `N = 4`

`XU` vector = `[x0; u0; x1; u1; x2; u2; x3; u3; x4]`

which is of size `((4 + 2) + (4 + 2) + (4 + 2) + (4 + 2) + (4)) = 28`
or equivalently `3Nn + 2n = 3(4)(2) + 2(2) = 24 + 4 = 28`, as desired.
"""
function initializeTraj(x0::Array{Float64, 1}, xN::Array{Float64, 1},
                        NSteps::Int64)
    if size(x0) != size(xN)
        error("x0 and xN must have the same size. " *
                "Currently: x0 = $x0 and xN = $xN")
    end
    if mod(size(x0, 1), 2) != 0
        error("x0 and xN must have an even number of variables. " *
              "Example, x0 = [sx0; sy0; vx0; vy0]. Currently x0 = $x0.")
    end

    # Calculate the linearly interpolated x values
    xSteps = LinRange(x0, xN, NSteps + 1)
    # Now we want to stack this as a vector in the form
    # [x0; u0; x1; ...; uN-1, xN]
    nDim = Int64(size(x0, 1) / 2)
    sizeXU = 3 * NSteps * nDim + 2 * nDim
    XUfull = zeros(sizeXU, 1)

    # println("Dimensions = $nDim")

    for k in 1:(NSteps + 1)
        iStart = 1 + (k - 1) * (3 * nDim)
        iEnd = iStart + (2 * nDim - 1)
        # println("Accessing ($iStart, $iEnd)")
        XUfull[iStart:iEnd] = xSteps[k]
    end

    return XUfull
end

function initializeTraj(x0::Array{Int64, 1}, xN::Array{Int64, 1},
                        NSteps::Int64)
    return initializeTraj(float.(x0), float.(xN), NSteps)
end

"""
initializeTraj(x0::Array{Float64, 1}, xN::Array{Float64, 1},
               u0::Array{Float64, 1}, uN::Array{Float64, 1},
               NSteps::Int64)

Initializes the full `XU` vector (see example below) using a linear
interpolation from the initial point `x0` to the final point `xN`.

Take `x0 = [s0; v0]` and `xN = [sN; vN`] to be of size (`2n×1`). The full `XU`
vector size is `(3Nn + 2n)×1`.

For example, take a 2D system: `x0 = [sx0; sy0; vx0; vy0]` and
                               `xN = [sxN; syN; vxN, vyN]`.

The controls are `uk = [uxk; uyk]`. Thus, for `N = 4`

`XU` vector = `[x0; u0; x1; u1; x2; u2; x3; u3; x4]`

which is of size `((4 + 2) + (4 + 2) + (4 + 2) + (4 + 2) + (4)) = 28`
or equivalently `3Nn + 2n = 3(4)(2) + 2(2) = 24 + 4 = 28`, as desired.
"""
function initializeTraj(x0::Array{Float64, 1}, xN::Array{Float64, 1},
                        u0::Array{Float64, 1}, uN::Array{Float64, 1},
                        NSteps::Int64)
    if size(x0) != size(xN)
        error("x0 and xN must have the same size. " *
                "Currently: x0 = $x0 and xN = $xN")
    end
    if mod(size(x0, 1), 2) != 0
        error("x0 and xN must have an even number of variables. " *
              "Example, x0 = [sx0; sy0; vx0; vy0]. Currently x0 = $x0.")
    end

    # Calculate the linearly interpolated x values
    xSteps = LinRange(x0, xN, NSteps + 1)
    uSteps = LinRange(u0, uN, NSteps)
    # Now we want to stack this as a vector in the form
    # [x0; u0; x1; ...; uN-1, xN]
    nDim = Int64(size(x0, 1) / 2)
    sizeXU = 3 * NSteps * nDim + 2 * nDim
    XUfull = zeros(sizeXU, 1)

    # println("Dimensions = $nDim")

    for k in 1:(NSteps + 1)
        iStart = 1 + (k - 1) * (3 * nDim)
        iEnd = iStart + (2 * nDim - 1)
        # println("Accessing ($iStart, $iEnd)")
        XUfull[iStart:iEnd] = xSteps[k]

        if k ≤ NSteps
            jStart = iEnd + 1
            jEnd = jStart + (nDim - 1)
            XUfull[jStart:jEnd] = uSteps[k]
        end
    end

    return XUfull
end
