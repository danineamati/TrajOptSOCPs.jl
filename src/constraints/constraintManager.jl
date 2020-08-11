# Constraint Manager holds a list of constaints and can run the necessary
# calculus

include("constraints.jl")

using SparseArrays

abstract type constraintManager end

"""
    constraintManager_Base

Contains a list of constraints and the dual variables (lambda) that correspond
with each constraint.

See also [`constraintManager_Dynamics`](@ref)
"""
struct constraintManager_Base <: constraintManager
    cList::Array{constraint, 1}
    lambdaList
end

#############################################################
###      Update the dual variable of the constraints      ###
#############################################################
@doc raw"""
    updateDual(cM::constraintManager, yCurr, penalty::Float64)

Updates (in place) each λ in the lambdaList stores in the
[`constraintManager`](@ref).

```math
λ ⟵ λ + ρ c(y)
```

for each constraint

"""
function updateDual!(cM::constraintManager, yCurr, penalty::Float64)
    for (i, c) in enumerate(cM.cList)
        # Obtain the dual variable matching this constraint
        lambda = cM.lambdaList[i]

        # Evaluate the current constraint
        if (typeof(c) == AL_simpleCone) | (typeof(c) == AL_Multiple_simpleCone)
            cVal = getNormToProjVals(c, yCurr, lambda)
        else
            cVal = getNormToProjVals(c, yCurr)
        end
        # cVal = getNormToProjVals(c, yCurr)

        lambda = lambda + penalty * cVal
    end
end

#############################################################
###            Evaluate a List of Constraints             ###
#############################################################

@doc raw"""
    evalConstraints(yCurr, cM::constraintManager, penalty::Float64)

Evaluates each constraint and sums the result.

Recall that the constraint terms are given as
```math
ϕ(x) = f(x) + ρ Σ ||cᵢ(x)||² + Σ λᵢ^{\top}cᵢ(x)
```

This function returns the following part:
```math
ρ Σ ||cᵢ(x)||² + Σ λᵢ^{\top}cᵢ(x)
```

See also: [`constraintManager`](@ref)
"""
function evalConstraints(cM::constraintManager, yCurr, penalty::Float64)

    total_eval = 0

    # println("Evaluating the constraints")

    for (i, c) in enumerate(cM.cList)
        # Obtain the dual variable matching this constraint
        lambda = cM.lambdaList[i]
        # Evaluate the current constraint
        if (typeof(c) == AL_simpleCone) | (typeof(c) == AL_Multiple_simpleCone)
            cVal = getNormToProjVals(c, yCurr, lambda)
        else
            cVal = getNormToProjVals(c, yCurr)
        end
        # cVal = getNormToProjVals(c, yCurr)
        # println("Lambda = $(size(lambda)) and cVal = $(size(cVal))")
        # println("So, ||c(x)|| = $(norm(cVal)) -> $(size(norm(cVal)))")
        # println("And λ'c = $(lambda' * cVal) -> $(size(lambda' * cVal))")
        # println()

        # Evaluate and add the full constraint term
        total_eval += penalty * (norm(cVal)^2) + (lambda' * cVal)[1]
    end

    return total_eval
end

"""
    getConstraintViolationList(cM::constraintManager, yList,
                               penalty::Float64 = 1.0)

Takes a list of trajectories and calculates the TOTAL constraint violation
at each trajectory.
"""
function getConstraintViolationList(cM::constraintManager, yList,
                                    penalty::Float64 = 1.0)
    constVio = zeros(size(yList, 1))

    for (ind, y) in enumerate(yList)
        constVio[ind] = evalConstraints(cM, y, penalty)
    end

    return constVio
end

#############################################################
###     Evaluate the Gradient of a List of Constraints    ###
#############################################################
@doc raw"""
    evalConstraints(yCurr, cM::constraintManager, penalty::Float64)

Evaluates the gradient of each constraint and sums the result.

Recall that the constraint terms are given as
```math
∇ϕ(x) = ∇f(x) + Σ J(c)'(ρ * cᵢ(x) + λᵢ)
```

This function returns the following part:
```math
Σ J(c)'(ρ * cᵢ(x) + λᵢ)
```

See also: [`constraintManager`](@ref)
"""
function evalGradConstraints(cM::constraintManager, yCurr, penalty::Float64)

    if size(yCurr, 1) == 1
        # If this is a one dimensional optimization problem, the gradient
        # is actually just a derivative. Hence the vector notation will cause
        # an error.
        total_grad = 0
    else
        # otherwise the vector notation is necessary
        total_grad = zeros(size(yCurr))
    end

    for (i, c) in enumerate(cM.cList)
        # Select the correct dual variables
        lambda = cM.lambdaList[i]
        # Evaluate the constraint
        if (typeof(c) == AL_simpleCone) | (typeof(c) == AL_Multiple_simpleCone)
            cVal = getNormToProjVals(c, yCurr, lambda)
        else
            cVal = getNormToProjVals(c, yCurr)
        end
        # cVal = getNormToProjVals(c, yCurr)
        # Evaluate the jacobian of the constraint
        cJacob = getGradC(c, yCurr)
        # println("cJacob:")
        # display(cJacob)
        # println("cVal:")
        # println(cVal)

        # Evaluate and add the full constraint term
        total_grad += cJacob' * (penalty * cVal + lambda)
    end

    return total_grad

end


#############################################################
###     Evaluate the Hessian of a List of Constraints     ###
#############################################################
@doc raw"""
    evalHessConstraints(yCurr, cM::constraintManager, penalty::Float64)

Evaluates the Hessian of each constraint and sums the result.

Recall that the constraint terms are given as
```math
H(ϕ(x)) = H(f(x)) + H(ρ Σ ||cᵢ(x)||² + Σ λᵢ^{\top}cᵢ(x))
```

This function returns the following part:
```math
H(ρ Σ ||cᵢ(x)||² + Σ λᵢ^{\top}cᵢ(x)) = Σ H(ρ ||cᵢ(x)||² + λᵢ^{\top}cᵢ(x))
```

See also: [`constraintManager`](@ref)
"""
function evalHessConstraints(cM::constraintManager, yCurr, rho = 1)

    total_hess = spzeros(size(yCurr, 1), size(yCurr, 1))

    for c in cM.cList
        total_hess += Symmetric(getHessC_ALTerm(c, yCurr, rho))
        #getHessC(c)#getHessC_ALTerm(c, yCurr, rho)
    end

    return total_hess

end
