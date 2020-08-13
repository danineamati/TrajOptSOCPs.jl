# SOCP Constraints for angle constraints

#=

These are constraints of the form
||[x1; x2]|| ≤ α x3
Since this is a trajectory, we also likely want to use this constraint
multiple times (e.g. at each step in the discretization)

This file handles these constraints such that it can be integrated into
the "constraintManager" framework.


Mathematics:
For the rocket soft-landing problem, we want to limit the angle that
engine can swing. We make the following assumptions:
1. The angle is usually to a small angle < 90°.
This holds because we don't want the engine gimbal to move too far from its
nominal axial position (i.e. engine exhause aligned with the rocket's body).
We also don't want the rocket control to go unstable due to a large torque
from a large angle.
2. The thrust vector is assumed to be _exactly_ in line with the engine
pointing vector.
This is a good approximation for a good rocket. There may be small variations,
but the engine is designed to direct the exhaust.
3. The rocket is oriented with the ground. (i.e. the long axial direction of
the rocket is perpendicular to the ground)
This is a good approximation near landing since the rocket landing gear must
downward facing to interface with the ground.

Then, we calculate the radial distance as ||ux|| in 2D and ||[ux; uy]|| in 3D.
The angle against the z component is
θ = arctan(||uxy|| / -uz)
The "-uz" is to have "-uz > 0" since the engine points downward. We want
θ < θ_max
Equivalently,
arctan(||uxy|| / -uz) < θ_max
||uxy|| < - tan(θ_max) uz

Call tan(θ_max) = α. Thus, we recover
||ux|| < - α uy in 2D
||[ux; uy]|| < - α uz in 3D

as desired.

=#


using SparseArrays, LinearAlgebra

include("projections.jl")


# -----------------------
# Simple Second-Order Cone Constraint for angle type constraints
# -----------------------
@doc raw"""
    AL_simpleAngleCone(alpha::Float64,
                       sIndicate::Diagonal{Int64,Array{Int64,1}},
                       tIndicate::Array{Int64,1})

Creates a simple Second Order Cone Constraint (e.g. for an SOCP) of the form
```math
||sIndicate * x|| ≤ alpha * tIndicate * x
```

So, in more standard, this is becomes:
```math
||sIndicate * x|| - alpha * tIndicate * x ≤ 0
```

For compactness, we will write this as
```math
||s|| - alpha t ≤ 0
```

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate the standard form.
- [`getNormToProjVals`](@ref) to evaluate the cone using projections.

"""
struct AL_simpleAngleCone <: constraint
    alpha::Float64
    sIndicate::Diagonal{Int64,Array{Int64,1}}
    tIndicate::Array{Int64,1}
end

get_s(r::AL_simpleAngleCone, x) = r.sIndicate * x
get_t(r::AL_simpleAngleCone, x) = r.alpha * r.tIndicate' * x


"""
    getRaw(r::AL_simpleAngleCone, x)

Evaluate the constraint without projections.
"""
getRaw(r::AL_simpleAngleCone, x) = norm(get_s(r, x)) - get_t(r, x)

"""
    satisfied(r::AL_simpleAngleCone, x)

Check constraint satisfaction. See also [`whichSatisfied`](@ref)
"""
satisfied(r::AL_simpleAngleCone, x) = (getRaw(r, x) ≤ 0)

"""
    whichSatisfied(r::AL_simpleAngleCone, xArr)

Check constraint satisfaction for each element in `xArr`.
"""
whichSatisfied(r::AL_simpleAngleCone, xArr) = [satisfied(r, x) for x in xArr]

"""
    getProjVecs(r::AL_simpleAngleCone, x, filled = true)

Gets the projection vectors for `x` onto the the second order cone defined
by `||sIndicate * x|| - alpha * tIndicate * x`.

See also [`AL_simpleAngleCone`](@ref)
"""
function getProjVecs(r::AL_simpleAngleCone, x, filled = true)
    return projSecondOrderCone(get_s(r, x), get_t(r, x), filled)
end

"""
    getNormToProjVals(r::AL_simpleAngleCone, x, λ = 0)

Get the constraint violation for the second order cone constraints.

The parameter `λ` describes the dual of the constraint and is used to determine
if the constraint is active. See [`coneActive`](@ref) for more information.
"""
function getNormToProjVals(r::AL_simpleAngleCone, x, λ = 0)
    s = get_s(r, x)
    t = get_t(r, x)
    active, filled = coneActive(s, t, λ)
    pv = getProjVecs(r, x, filled)

    return [norm([s; t] - pv)]
end


"""
    getGradC(r::AL_simpleAngleCone, t)

Calculates the gradient of a constraint.

For `r::AL_simpleAngleCone`, this is the `∇c(x) = [s / ||s||; - alpha]`

There is a singularity at the tip.
"""
function getGradC(r::AL_simpleAngleCone, x)
    s = get_s(r, x)

    # ds will have 0 in the disregarded spot for dt
    # Ex: s1 = Diagonal([1; 0; 1]) * [4; 5; 6] -> [4; 0; 6]
    #     ds1 = [0.5547001962252291;  0.0;  0.8320502943378437]
    ds = s / norm(s)

    # dt will have 0 in the disregarded spots for ds
    # Ex: t1 = [0; 1; 0]' * [4; 5; 6]
    #     dt1 = [0.0;  -0.2;  -0.0]
    dt = - r.alpha * r.tIndicate

    # Adding it together gets the actual grad.
    return ds + dt
end

"""
    getHessC(r::AL_simpleAngleCone, t)

Calculate the hessian of a constraint.

For `r::AL_simpleAngleCone`, `H(c(x)) = 0`
"""
function getHessC(r::AL_simpleAngleCone, x)
    return spzeros(size(x, 1), size(x, 1))
end

"""
    getHessC_ALTerm(r::AL_simpleAngleCone, t, rho = 1)

This is the hessian of the *constraint term* as it appears in the
augmented lagrangian.

For `r::AL_simpleAngleCone`,
`H(ρ c(x)'c(x) + λ c(x)) = ρ [s * s' / norm(s)^2, -alpha s/norm(s);
-alpha s/norm(s), alpha^2`
"""
function getHessC_ALTerm(r::AL_simpleAngleCone, x, rho = 1)
    #=
    The augmented lagrangian constraint term is of the form:
    ρ c(x)'c(x) + λ c(x)

    where c(x) = x / ||x||.

    The Hessian matrix is then
    Htot = ρ (c.H + J.J + H.c) + λ H

    But H = 0 for Affine equalities, so when the constraint is not satisfied,
    Htot = ρ (J.J) =
    ρ * [s * s' / norm(s)^2     -alpha s/norm(s)]
        [-alpha s'/norm(s)       alpha^2]
    When the constraint is satisfied, it is just zero.

    However, we want this to work which ever order we have "s" and "t". And,
    more precisely, we have "s" as the same length as "x", with zeros in the
    gaps for t. Hence, we use four indicators:

    s * s' / norm(s)^2 -> sI2 * sI2'
    -alpha s/norm(s)   -> sI2 * tI'
    -alpha s'/norm(s)  -> (sI2 * tI')'
    alpha^2            -> tI * tI'

    where tI is r.tIndicate and sI2 is diag(sI)
    =#

    if satisfied(r, x)
        # Returns a sparse zero matrix of the right dimensions
        return getHessC(r, x)
    end

    # Get s. t is not used
    s = get_s(r, x)
    ns = norm(s)

    # Get the indicators
    sI2 = diag(r.sIndicate)
    tI = r.tIndicate

    a = (s * s' / (ns^2)) .* (sI2 * sI2')
    b = (- r.alpha * s / ns) .* (sI2 * tI')
    bT = (- r.alpha * s' / ns) .* (sI2 * tI')'
    c = (r.alpha^2) * (tI * tI')

    return rho * Symmetric(a + b + bT + c)
end

# -----------------------
# Simple Second-Order Cone Constraints for Multiple angle type constraints
# -----------------------
@doc raw"""
    AL_Multiple_AngleCone(
                        sc::AL_simpleAngleCone,
                        nDim::Int64,
                        indicatorList::Tuple{Array{Int64,1},Array{Int64,1}})

Handles Multiple Second Order Cone Constraint (SOCP) of the form
```math
||s_k|| ≤ alpha t_k \quad \text{for} \ k ∈ N
```
where T is a list of primal decision vectors.

So, in more standard form, this is becomes:
```math
||s_k|| - alpha t_k ≤ 0 \quad \text{for} \ k ∈ N
```

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate each cone without projections.
- [`getNormToProjVals`](@ref) to evaluate each cone using projections.

"""
struct AL_Multiple_AngleCone <: constraint
    sc::AL_simpleAngleCone
    nDim::Int64
    #Type of findnz(sparseArray)[1] or findall(!iszero, sparseArray)
    indicatorList::Array{Int64,1}
end

function makeAL_Multiple_AngleCone(alpha::Float64,
                sIndicate::Diagonal{Int64,Array{Int64,1}},
                tIndicate::Array{Int64,1}, nDim::Int64,
                indicatorList::Array{Int64,1})
    return AL_Multiple_AngleCone(
                    AL_simpleAngleCone(alpha, sIndicate, tIndicate),
                    nDim, indicatorList)
end

"""
    parseRelevantSteps(r::AL_Multiple_AngleCone, x)

Run through the indicator list and use the dimensions to extract the
relevant vectors from x
"""
function parseRelevantSteps(r::AL_Multiple_AngleCone, x)
    # This array will store the relevant parts of the x vector.
    steps = []

    for ind in r.indicatorList
        # Find the relevant parts and push them as a group
        push!(steps, x[ind:(ind + r.nDim - 1)])
    end

    return steps
end

"""
    getRaw(r::AL_Multiple_AngleCone, x)

Evaluate each of the simple cone constraints without projection.
"""
function getRaw(r::AL_Multiple_AngleCone, x)
    # Use the parsed steps and calculate the raw value for each set.
    steps = parseRelevantSteps(r, x)
    return [getRaw(r.sc, s) for s in steps]
end


"""
    getProjVecs(r::AL_Multiple_AngleCone, steps, filled)

Get the projection defined by [`AL_simpleAngleCone`](@ref) for each of the
constraints.

See also [`AL_Multiple_AngleCone`](@ref)
"""
function getProjVecs(r::AL_Multiple_AngleCone, steps, filled)
    # Then loop through the parsed steps to get the projection.
    pv = []
    for (i, s) in enumerate(steps)
        # Filled is needed to determine which projection and may be different
        # for each of the different SOCP constraints in question.
        push!(pv, getProjVecs(r.sc, s, filled[i]))
    end

    return pv
end

"""
    getNormToProjVals(r::AL_Multiple_AngleCone, x, λ)

Get the constraint violation for multiple second order cone constraints. The
length of the vector `λ` must match the number of steps.

See also [`AL_Multiple_AngleCone`](@ref) and [`coneActive`](@ref) for more
information.
"""
function getNormToProjVals(r::AL_Multiple_AngleCone, x, λ)
    # Use the parsed steps to separate our each constraint
    steps = parseRelevantSteps(r, x)

    # Determine which of the constraints (if any) require an unfilled cone/
    fillVec = Bool[]
    for (ind, st) in enumerate(steps)
        s = get_s(r.sc, st)
        t = get_t(r.sc, st)
        push!(fillVec, coneActive(s, t, λ[ind])[2])
    end

    # Using the filled information and the steps in question, get the
    # projection vectors
    pvs = getProjVecs(r, steps, fillVec)

    # Use the projection vectors to calculate the constraints.
    constVio = []
    for (st, p) in zip(steps, pvs)
        s = get_s(r.sc, st)
        t = get_t(r.sc, st)
        push!(constVio, norm([s; t] - p))
    end

    # If there are NaN values (this happens if you take the norm of a zero
    # vector that is large), use the snippet of code below.
    return replace(constVio, NaN=>0.0)
end

function getNormToProjVals(r::AL_Multiple_AngleCone, x)
    lambdaVec = spzeros(size(r.indicatorList, 1))
    return getNormToProjVals(r, x, lambdaVec)
end
