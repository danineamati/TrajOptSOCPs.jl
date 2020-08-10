# Only the affine constraints

using SparseArrays

include("..\\other_utils\\feasibleCheck.jl")
include("projections.jl")

# ---------------------
# Equality constraints
# ---------------------
"""
    AL_AffineEquality(A, b)

Creates an Affine Equality Constraint of the form `Ax = b`.

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate `Ax - b`
- [`getNormToProjVals`](@ref) to evaluate the projection to `Ax = b`
"""
struct AL_AffineEquality <: constraint
    A
    b
end

# Check that a constraint is satisfied (Ax - b = 0)
"""
    satisfied(r::AL_AffineEquality, x) = (r.A * x == r.b)

Check constraint satisfaction. See also [`whichSatisfied`](@ref)
"""
satisfied(r::AL_AffineEquality, x) = (r.A * x == r.b)

"""
    whichSatisfied(r::AL_AffineEquality, x) = (r.A * x .== r.b)

Check constraint satisfaction. See also [`satisfied`](@ref)
"""
whichSatisfied(r::AL_AffineEquality, x) = (r.A * x .== r.b)

# Evaluate the constraint.
"""
    getRaw(r::AL_AffineEquality, x) = r.A * x - r.b

Evaluate the constraint without projection
"""
getRaw(r::AL_AffineEquality, x) = r.A * x - r.b

"""
    getProjVecs(r::AL_AffineEquality, x)

For every row, get the projection onto the line `a'x = b`.

See [`projAffineEq`](@ref) for more information.
"""
function getProjVecs(r::AL_AffineEquality, x)
    #=
    We want to project each constraint
    We write Ax = b
    But this is equivalent to a1'x = b1, a2'x = b2, ..., am'x = bm
    Where ak is the row k in A
    =#
    projVecs = []
    for ind in 1:size(r.b, 1)
        aI = r.A[ind, :]
        bI = r.b[ind]
        push!(projVecs, projAffineEq(aI, bI, x))
    end

    return projVecs
end

"""
    getNormToProjVals(r::AL_AffineEquality, x)

Get the row by row constraint violation based on projections. This is better
than [`getRaw`](@ref) since it doesn't depend on the magnitude of the elements
in `r.A` and `r.b`.

See also [`AL_AffineEquality`](@ref), [`getProjVecs`](@ref),
[`projAffineEq`](@ref)
"""
function getNormToProjVals(r::AL_AffineEquality, x)
    #=
    We want to get the projected vector and calculate the distance between the
    original point and the constraint.
    =#
    projVecs = getProjVecs(r, x)
    pvsDiff = [pv - x for pv in projVecs]

    # NaN removal due to floating point error near 0.0
    return replace(norm.(pvsDiff, 2), NaN=>0.0)
end

# Lastly, we do some calculus
"""
    getGradC(r::AL_AffineEquality, x)

Calculates the gradient of a constraint.

For `r::AL_AffineEquality`, `∇c(x) = A`

The variable "x" is passed to the header to maintain the same function header
but it is not used in the function.
"""
function getGradC(r::AL_AffineEquality, x)
    #=
    We want the gradient of the constraints. This is piecewise in the
    inequality case. But it is a single function in the equality case
    =#
    return r.A
end

"""
    getHessC(r::AL_AffineEquality)

Calculates the hessian of a constraint.

For `r::AL_AffineEquality`, `H(c(x)) = 0`
"""
function getHessC(r::AL_AffineEquality)
    #=
    We want the hessian of the constraints. This is just zero for affine
    constraints, but we need to match the dimensions.
    =#
    return spzeros(size(r.A, 2), size(r.A, 2)) # nxn
end

"""
    getHessC_ALTerm(r::AL_AffineEquality, x, rho = 1)

Calculates the hessian of the augmented lagrangian constraint term.

For `r::AL_AffineEquality`, `H(ρ c(x)'c(x) + λ c(x)) = ρ A'A`

The variable "x" is passed to the header to maintain the same function header
but it is not used in the function.
"""
function getHessC_ALTerm(r::AL_AffineEquality, x, rho = 1)
    #=
    The augmented lagrangian constraint term is of the form:
    ρ c(x)'c(x) + λ c(x)

    where c(x) = Ax - b.

    The Hessian matrix is then
    Htot = ρ (c.H + J.J + H.c) + λ H

    But H = 0 for Affine equalities, so
    Htot = ρ (J.J) = ρ A'A
    =#

    jacob = getGradC(r, x)
    return rho * jacob'jacob
end


# -----------------------
# Inequality constraints
# -----------------------

"""
    AL_AffineInequality(A, b)

Creates an Affine Inequality Constraint of the form `Ax ≤ b`.

To check constraint satisfaction, use:
[`satisified`](@ref) and [`whichSatisfied`](@ref)

To evaluate the constraint, use:
- [`getRaw`](@ref) to evaluate `Ax - b`
- [`getNormToProjVals`](@ref) to evaluate the projection to `Ax ≤ b`
"""
struct AL_AffineInequality <: constraint
    A
    b
end

# Check that a constraint is satisfied (Ax - b ≤ 0)
"""
    satisfied(r::AL_AffineInequality, x)

Check constraint satisfaction. See also [`whichSatisfied`](@ref)
"""
satisfied(r::AL_AffineInequality, x) = isFeasiblePolyHedron(r.A, r.b, x)

"""
    whichSatisfied(r::AL_AffineInequality, x)

Check constraint satisfaction. See also [`satisfied`](@ref)
"""
whichSatisfied(r::AL_AffineInequality, x) = (r.A * x .≤ r.b)

# Evaluate the constraint.
"""
    getRaw(r::AL_AffineInequality, x) = r.A * x - r.b

Evaluate the constraint without projection
"""
getRaw(r::AL_AffineInequality, x) = r.A * x - r.b

"""
    getProjVecs(r::AL_AffineInequality, x)

For every row, get the projection onto the line `a'x = b`.

See [`projAffineEq`](@ref) for more information.
"""
function getProjVecs(r::AL_AffineInequality, x)
    #=
    We want to project each constraint
    We write Ax = b
    But this is equivalent to a1'x = b1, a2'x = b2, ..., am'x = bm
    Where ak is the row k in A
    =#
    projVecs = []
    for ind in 1:size(r.b, 1)
        aI = r.A[ind, :]
        bI = r.b[ind]
        push!(projVecs, projAffineIneq(aI, bI, x))
    end

    return projVecs
end

"""
    getNormToProjVals(r::AL_AffineInequality, x)

Get the row by row constraint violation based on projections. This is better
than [`getRaw`](@ref) since it doesn't depend on the magnitude of the elements
in `r.A` and `r.b`. It also automatically accounts for the inequality.

See also [`AL_AffineInequality`](@ref), [`getProjVecs`](@ref),
[`projAffineEq`](@ref)
"""
function getNormToProjVals(r::AL_AffineInequality, x)
    #=
    We want to get the projected vector and calculate the distance between the
    original point and the constraint.
    =#
    projVecs = getProjVecs(r, x)
    pvsDiff = [pv - x for pv in projVecs]
    return norm.(pvsDiff, 2)
end

# Lastly, we do some calculus
"""
    getGradC(r::AL_AffineInequality, x)

Calculates the gradient of a constraint.

For `r::AL_AffineInequality`, `∇c(x) = A` for rows where the constraint is NOT
satisfied.

Unlike the [`AL_AffineEquality`](@ref) case, the variable "x" is passed to
is needed to determine if the constraint is satisfied.
"""
function getGradC(r::AL_AffineInequality, x)
    #=
    We want the gradient of the constraints. This is piecewise in the
    inequality case. But it is a single function in the equality case
    =#
    APost = zeros(size(r.A))
    rowPassed = whichSatisfied(r, x)
    for row in 1:size(r.A, 1)
        if !rowPassed[row]
            # Constraint is not met
            APost[row, :] = r.A[row, :]
        end
    end
    return APost
end

"""
    getHessC(r::AL_AffineEquality)

Calculates the hessian of a constraint.

For `r::AL_AffineInequality`, `H(c(x)) = 0`
"""
function getHessC(r::AL_AffineInequality)
    #=
    We want the hessian of the constraints. This is just zero for affine
    constraints, but we need to match the dimensions.
    =#
    return spzeros(size(r.A, 2), size(r.A, 2)) # nxn
end

"""
    getHessC_ALTerm(r::AL_AffineInequality, x, rho = 1)

Calculates the hessian of the augmented lagrangian constraint term.

For `r::AL_AffineEquality`, `H(ρ c(x)'c(x) + λ c(x)) = ρ A+'A+`

The variable "x" is passed to the header to maintain the same function header
but it is not used in the function.
"""
function getHessC_ALTerm(r::AL_AffineInequality, x, rho = 1)
    #=
    The augmented lagrangian constraint term is of the form:
    ρ c(x)'c(x) + λ c(x)

    where c(x) = Ax - b.

    The Hessian matrix is then
    Htot = ρ (c.H + J.J + H.c) + λ H

    But H = 0 for Affine equalities, so
    Htot = ρ (J.J) = ρ A'A
    =#

    jacob = getGradC(r, x)
    return sparse(rho * jacob'jacob)
end
