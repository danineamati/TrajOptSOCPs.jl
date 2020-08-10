# Only NEW SOCP Constraints in comparison to the original SOCP solver

using SparseArrays

include("..\\other_utils\\feasibleCheck.jl")
include("projections.jl")


# -----------------------
# Simple Second-Order Cone Constraints
# -----------------------
@doc raw"""
    AL_coneSlack(tMax)

Creates a _simple_ Second Order Cone Constraint (SOCP) of the form
```math
||t|| ≤ tMax
```

So, in more standard form, this is becomes:
```math
||t|| - tMax ≤ 0
```

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate `||t|| - tMax`.
- [`getNormToProjVals`](@ref) to evaluate the projection to the affine
equalities and (more critically) the second order cone.

"""
struct AL_simpleCone <: constraint
    tMax
end

"""
    getRaw(r::AL_simpleCone, t)

Evaluate the constraint without projection.
"""
getRaw(r::AL_simpleCone, t) = norm(t) - r.tMax

"""
    satisfied(r::AL_simpleCone, t)

Check constraint satisfaction. See also [`whichSatisfied`](@ref)
"""
satisfied(r::AL_simpleCone, t) = (getRaw(r, t) ≤ 0)

"""
    whichSatisfied(r::AL_simpleCone, tArr)

Check constraint satisfaction for each element in `tArr`.
"""
whichSatisfied(r::AL_simpleCone, tArr) = [satisfied(r, t) for t in tArr]

"""
    getProjVecs(r::AL_simpleCone, t, filled = true)

Gets the projection vectors for `t` onto the the second order cone defined
by `||t|| - tMax`.

See also [`AL_simpleCone`](@ref)
"""
function getProjVecs(r::AL_simpleCone, t, filled = true)
    return projSecondOrderCone(t, r.tMax, filled)
end

"""
    getNormToProjVals(r::AL_simpleCone, t, λ = 0)

Get the constraint violation for the second order cone constraints.

The parameter `λ` describes the dual of `||t|| - tMax` and is used to determine
if the constraint is active. See [`coneActive`](@ref) for more information.
"""
function getNormToProjVals(r::AL_simpleCone, t, λ = 0)

    active, filled = coneActive(t, r.tMax, λ)
    pv = getProjVecs(r, t, filled)

    return [norm([t; r.tMax] - pv)]
end

"""
    getGradC(r::AL_simpleCone, t)

Calculates the gradient of a constraint.

For `r::AL_simpleCone`, this is the `∇c(x) = x / ||x||`

There is a singularity at the tip.
"""
function getGradC(r::AL_simpleCone, t)
    return t / norm(t)
end

"""
    getHessC(r::AL_simpleCone, t)

Calculate the hessian of a constraint.

For `r::AL_simpleCone`, `H(c(x)) = 0`
"""
function getHessC(r::AL_simpleCone, t)
    return spzeros(size(t, 1), size(t, 1))
end

"""
    getHessC_ALTerm(r::AL_simpleCone, t, rho = 1)

This is the hessian of the *constraint term* as it appears in the
augmented lagrangian.

For `r::AL_simpleCone`, `H(ρ c(x)'c(x) + λ c(x)) = ρ s * s' / norm(s)^2`
"""
function getHessC_ALTerm(r::AL_simpleCone, t, rho = 1)
    #=
    The augmented lagrangian constraint term is of the form:
    ρ c(x)'c(x) + λ c(x)

    where c(x) = x / ||x||.

    The Hessian matrix is then
    Htot = ρ (c.H + J.J + H.c) + λ H

    But H = 0 for Affine equalities, so when the constraint is not satisfied,
    Htot = ρ (J.J) = ρ s * s' / norm(s)^2
    Wehen the constraint is satisfied, it is just zero.
    =#

    if satisfied(r, t)
        # Returns a sparse zero matrix of the right dimensions
        return getHessC(r, t)
    end

    return rho * (t * t') / (norm(t)^2)
end


# -----------------------
# Simple Second-Order Cone Constraints for Multiple
# -----------------------
@doc raw"""
    AL_Multiple_simpleCone(
                        sc::AL_simpleCone,
                        nDim::Int64,
                        indicatorList::Tuple{Array{Int64,1},Array{Int64,1}})

Handles Multiple Second Order Cone Constraint (SOCP) of the form
```math
||t_k|| ≤ t_{Max} \quad \text{for} \ t_k ∈ T
```
where T is a list of primal decision vectors.

So, in more standard form, this is becomes:
```math
||t_k|| - t_{Max} ≤ 0 \quad \text{for} \ t_k ∈ T
```

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate `||t|| - tMax`.
- [`getNormToProjVals`](@ref) to evaluate the projection to the affine
equalities and (more critically) the second order cone.

"""
struct AL_Multiple_simpleCone <: constraint
    sc::AL_simpleCone
    nDim::Int64
    #Type of findnz(sparseArray)[1] or findall(!iszero, sparseArray)
    indicatorList::Array{Int64,1}
end

function makeAL_Multiple_simpleCone(tMax::Float64, nDim::Int64,
                                    indicatorList::Array{Int64,1})
    return AL_Multiple_simpleCone(AL_simpleCone(tMax), nDim, indicatorList)
end



"""
    parseRelevantSteps(r::AL_Multiple_simpleCone, x)

Run through the indicator list and use the dimensions to extract the
relevant vectors from x
"""
function parseRelevantSteps(r::AL_Multiple_simpleCone, x)
    # This array will store the relevant parts of the x vector.
    steps = []

    for ind in r.indicatorList
        # Find the relevant parts and push them as a group
        push!(steps, x[ind:(ind + r.nDim - 1)])
    end

    return steps
end

"""
    getRaw(r::AL_Multiple_simpleCone, x)

Evaluate each of the simple cone constraints without projection.
"""
function getRaw(r::AL_Multiple_simpleCone, x)
    # Use the parsed steps and calculate the raw value for each set.
    steps = parseRelevantSteps(r, x)
    return [getRaw(r.sc, s) for s in steps]
end

"""
    getProjVecs(r::AL_Multiple_simpleCone, steps, filled)

Get the projection defined by `||t|| - tMax` for each of the constraints.

See also [`AL_Multiple_simpleCone`](@ref)
"""
function getProjVecs(r::AL_Multiple_simpleCone, steps, filled)
    # Then loop through the parsed steps to get the projection.
    pv = []
    for (i, s) in enumerate(steps)
        # Filled is needed to determine which projection and may be different
        # for each of the different SOCP constraints in question.
        push!(pv, getProjVecs(r.sc, s, filled[i]))
    end

    return pv
end

function getProjVecs(r::AL_Multiple_simpleCone, x::Array{Float64, 1}, filled)
    steps = parseRelevantSteps(r, x)
    return getProjVecs(r, steps, filled)
end

function getProjVecs(r::AL_Multiple_simpleCone, x::Array{Int64, 1}, filled)
    steps = parseRelevantSteps(r, x)
    return getProjVecs(r, steps, filled)
end

"""
    getNormToProjVals(r::AL_Multiple_simpleCone, x, λ)

Get the constraint violation for multiple second order cone constraints. The
length of the vector `λ` must match the number of steps.

See also [`AL_Multiple_simpleCone`](@ref) and [`coneActive`](@ref) for more
information.
"""
function getNormToProjVals(r::AL_Multiple_simpleCone, x, λ)
    # Use the parsed steps to separate our each constraint
    steps = parseRelevantSteps(r, x)

    # Determine which of the constraints (if any) require an unfilled cone/
    fillVec = Bool[]
    for (ind, s) in enumerate(steps)
        push!(fillVec, coneActive(s, r.sc.tMax, λ[ind])[2])
    end

    # Using the filled information and the steps in question, get the
    # projection vectors
    pvs = getProjVecs(r, steps, fillVec)

    # Use the projection vectors to calculate the constraints.
    constVio = []
    for (s, p) in zip(steps, pvs)
        push!(constVio, norm([s; r.sc.tMax] - p))
    end

    # If there are NaN values (this happens if you take the norm of a zero
    # vector that is large), use the snippet of code below.
    return replace(constVio, NaN=>0.0)
end

function getNormToProjVals(r::AL_Multiple_simpleCone, x)
    lambdaVec = spzeros(size(r.indicatorList, 1))
    return getNormToProjVals(r::AL_Multiple_simpleCone, x, lambdaVec)
end

"""
    getGradC(r::AL_Multiple_simpleCone, x)

Parses the relevant `x` vector for each cone constraint. It calculates the
gradient of each and then stiches it into a jacobian of the same dimension
as `x` and number of constraints with the zeros in relevant indicies.

(Result is a sparse array)

See also [`AL_Multiple_simpleCone`](@ref) and [`AL_simpleCone`](@ref)
"""
function getGradC(r::AL_Multiple_simpleCone, x)
    # Use the parsed steps to separate our each constraint
    steps = parseRelevantSteps(r, x)

    # Use the AL_simpleCone to get the gradient for each constraint
    gradList = [getGradC(r.sc, s) for s in steps]

    # Stitch it back together.
    gradStitched = spzeros(size(r.indicatorList, 1), size(x, 1))

    for (ind, pos) in enumerate(r.indicatorList)
        posEnd = pos + r.nDim - 1
        gradStitched[ind, pos:posEnd] = gradList[ind]
    end

    return gradStitched
end

"""
    getHessC(r::AL_Multiple_simpleCone, x)

Calculate the hessian of a constraint.

For `r::AL_Multiple_simpleCone`, `H(c(x)) = 0`

See also [`AL_Multiple_simpleCone`](@ref) and [`AL_simpleCone`](@ref)
"""
function getHessC(r::AL_Multiple_simpleCone, x)
    return spzeros(size(x, 1), size(x, 1))
end

"""
    getHessC_ALTerm(r::AL_Multiple_simpleCone, x, rho = 1)

This is the hessian of the *constraint term* as it appears in the
augmented lagrangian.

For `r::AL_Multiple_simpleCone`,
`H(ρ c(x)'c(x) + λ c(x)) = ρ s * s' / norm(s)^2` for each constraint.

Parses the relevant `x` vector for each cone constraint. It calculates the
hessian of each and then stiches it into a vector of the same dimension
as `x`×`x` with the zeros in relevant indicies.

(Result is a sparse array)

See also [`AL_Multiple_simpleCone`](@ref) and [`AL_simpleCone`](@ref)
"""
function getHessC_ALTerm(r::AL_Multiple_simpleCone, x, rho = 1)
    #=
    The augmented lagrangian constraint term is of the form:
    ρ c(x)'c(x) + λ c(x)

    where c(x) = x / ||x||.

    The Hessian matrix is then
    Htot = ρ (c.H + J.J + H.c) + λ H

    But H = 0 for Affine equalities, so
    Htot = ρ (J.J) = ρ s * s' / norm(s)^2

    For each constraint!
    =#

    # Use the parsed steps to separate our each constraint
    steps = parseRelevantSteps(r, x)

    # Use the AL_simpleCone to get the total hessian for each constraint
    hessList = [Symmetric(getHessC_ALTerm(r.sc, s)) for s in steps]

    # Stitch it back together.
    hessStitched = spzeros(size(x, 1), size(x, 1))

    for (ind, pos) in enumerate(r.indicatorList)
        posEnd = pos + r.nDim - 1
        hessStitched[pos:posEnd, pos:posEnd] = hessList[ind]
    end

    return hessStitched

end
