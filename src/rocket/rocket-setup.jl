# Contains the rocket set up and structs

#=
Now we want to set up the rocket dynamics. We have a linear system:
xk+1 = Axk + Buk + G

where xk = [sk; vk] (the position and velocity, respectively). If sk is n×1,
A is 2n×2n, B is 2n×n, and G is a 2n×1 vector.

Thus, we get the system:
[A B -I] * [xk; uk; xk+1] = G
which we can stack as

[A B -I 0 0 0 0 0 0]
[0 0 A B -I 0 0 0 0]
[0 0 0 0 A B -I 0 0]
[0 0 0 0 0 0 A B -I]
=#

using SparseArrays, LinearAlgebra

# include("trajectory-setup.jl")


"""
    rocket_simple(mass, isp, grav, deltaTime)

Barebones "spherical" rocket. This is to say, this struct only holds the mass
and the isp (specific impulse) of the rocket.
"""
struct rocket_simple
    mass
    isp
    grav
    deltaTime
end


@doc raw"""
    rocketDynamics(r::rocket_simple, nDim::Int64)

We have a linear dynamics system, that when continuous is
```math
\frac{dx}{dt} = Ax + Bu + g
```
If we have `y = [x; u; g]`, then,
```math
\frac{dy}{dt} = [A B I; 0 0 0; 0 0 0] y = W y
```
which has solution
```math
y(t) = e^{W t} y0
```
where W is square. In a discrete system,
```math
yk+1 = e^{W Δt} yk = E yk
```

Note that due to all the zeros in W, we get:
```math
xk+1 = E11 xk + E12 uk + E13 g
```
```math
uk+1 = E21 xk + E22 uk + E23 g = uk
```
```math
g = E31 xk + E32 uk + E33 g = g
```

Thus, once we calculate E, we extract E11, E12, and E13 to a system of the
form:
```math
x_{k+1} = Aₑ x_k + Bₑ u_k + Gₑ
```
We compute and return Aₑ, Bₑ, and Gₑ for a [`rocket_simple`](@ref).

For a simple rocket, we have

```math
A = [0 \ I; 0 \ 0] \quad B = [0; \frac{1}{m} I] \quad G = [0; -g]
```

We then take the exponent and return a resulting sparse matrices

"""
function rocketDynamics(r::rocket_simple, nDim::Int64)

    # Start off with a dense representation bc it is needed for the
    # exponentiation.
    z0 = zeros(nDim, nDim)
    aMat = [z0 I; z0 z0]
    bMat = [z0; -(1/r.mass) * I]
    gMat = [z0; I]
    # gVec = [spzeros(nDim, 1); r.grav]

    # Top Row
    abg = [aMat bMat gMat]
    abg_size = size(abg)

    # We need to pad the remaining rows to make a square matrix
    remainingRows = abg_size[2] - abg_size[1]
    abgTot = [abg; zeros(remainingRows, abg_size[2])]

    # Take the exponent e^(W Δt)
    eMat = exp(abgTot * r.deltaTime)

    # Extract the effective matrices
    rowsA = size(aMat, 1)
    colsA = size(aMat, 2)
    colsB = size(bMat, 2)

    aEffMat = sparse(eMat[1:rowsA, 1:colsA])
    bEffMat = sparse(eMat[1:rowsA, colsA + 1:colsA + colsB])
    gEffMat = sparse(eMat[1:rowsA, colsA + colsB + 1:end])

    return aEffMat, bEffMat, gEffMat
end

@doc raw"""
    rocketDynamicsStack(r::rocket_simple, nDim::Int64, NSteps::Int64)

We want to take the per discretization dynamics and write it as a full matrix.

Specifically, we get the system:
`[A B -I] * [xk; uk; xk+1] = G`
which we can stack as

```math
[A \ B \ -I \ 0 \ 0 \ 0 \ 0 \ 0 \ 0]
```

```math
[0 \ 0 \ A \ B \ -I \ 0 \ 0 \ 0 \ 0]
```

```math
[0 \ 0 \ 0 \ 0 \ A \ B \ -I \ 0 \ 0]
```

```math
[0 \ 0 \ 0 \ 0 \ 0 \ 0 \ A \ B \ -I]
```

We return the above matrix and the stacked `G` matrix.

See: [`rocketDynamicsFull`](@ref) and [`rocketDynamics`](@ref)
"""
function rocketDynamicsStack(r::rocket_simple, nDim::Int64, NSteps::Int64)
    # First calculate the A and B matrices. The C vector is not needed in this
    # case.
    Ak, Bk, Gk = rocketDynamics(r::rocket_simple, nDim::Int64)
    ABI_unit = [Ak Bk -I]
    colLen = size(ABI_unit, 2)

    # At each step we go down 2 * n rows and move right 3 * n columns.
    sizeCols = 3 * NSteps * nDim + 2 * nDim
    sizeRows = 2 * nDim * NSteps

    # Initialize and empty sparse A matrix and G vector
    AStacked = spzeros(sizeRows, sizeCols)
    GStacked = spzeros(sizeRows, 1)

    rStart = 1
    cStart = 1
    for k in 1:NSteps
        rEnd = rStart + (2 * nDim - 1)
        cEnd = cStart + colLen - 1

        # First, we calculate that stacked A Matrix
        # println("Accessing $rStart:$rEnd and $cStart:$cEnd")
        AStacked[rStart:rEnd, cStart:cEnd] = ABI_unit

        # Second, we calculate the stacked G vector
        GStacked[rStart:rEnd] = -Gk * r.grav

        # Prepare for next iteration.
        rStart += 2 * nDim
        cStart += 3 * nDim
    end

    return AStacked, GStacked
end


"""
    rocketDynamicsFull(r::rocket_simple, x0::Array{Float64, 1},
                             xN::Array{Float64, 1}, NSteps::Int64)

We combine the stack at each `xk` with the initial and final conditions. Namely,

`[I0; AStack; IN] * XU = [x0; GStack; xN]`

where `XU` is the stacked trajectory vector (x states and u controls). We
return `[I0; AStack; IN]` and `[x0; GStack; xN]`

See [`rocketDynamicsStack`](@ref), [`rocketDynamics`](@ref),
[`rocket_simple`](@ref)
"""
function rocketDynamicsFull(r::rocket_simple, x0::Array{Float64, 1},
                                 xN::Array{Float64, 1}, NSteps::Int64)
     if size(x0) != size(xN)
         error("x0 and xN must have the same size. " *
                 "Currently: x0 = $x0 and xN = $xN")
     end
     if mod(size(x0, 1), 2) != 0
         error("x0 and xN must have an even number of variables. " *
               "Example, x0 = [sx0; sy0; vx0; vy0]. Currently x0 = $x0.")
     end
    # The number of columns must be the same for all three. (Stacked vertically)
    nDim = Int64(size(x0, 1) / 2)
    sizeCols = 3 * NSteps * nDim + 2 * nDim
    sizeRowsI = 2 * nDim

    # Use the above function to determine the stacked matrix in the middle.
    AStacked, GStacked = rocketDynamicsStack(r, nDim, NSteps)
    I0 = [I spzeros(sizeRowsI, sizeCols - sizeRowsI)]
    IN = [spzeros(sizeRowsI, sizeCols - sizeRowsI) I]

    AFull = [I0; AStacked; IN]
    BFull = [x0; GStacked; xN]
    return AFull, BFull
end

function rocketDynamicsFull(r::rocket_simple, x0::Array{Int64, 1},
                                 xN::Array{Int64, 1}, NSteps::Int64)
    return rocketDynamicsFull(r, float.(x0), float.(xN), NSteps)
end
