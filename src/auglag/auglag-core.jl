# We have a similar QP set-up as before

# Simply include this file

# QP Set-Up
# By definition, the QP has the following Properties
# minimize_x (1/2) xT Q x + cT x
# subject to Ax ≤ B
#
# where (all are Real)
# Q is an nxn Symmetric Matrix (that is positive definite for convex)
# x is an nx1 vector
# c is an nx1 vector
# A is an mxn Matrix
# b is an mx1 vector

using LinearAlgebra, SparseArrays

include("..\\objective\\QP_Linear_objectives.jl")
include("..\\constraints\\constraintManager.jl")
include("..\\other_utils\\utils.jl")


# --------------------------
# Lagrangian
# φ(y) = f(x) + (ρ/2) ||c(y)||_2^2 + λ c(y)
#      = f(x) + (ρ/2) c(y)'c(y)    + λ c(y)
# --------------------------
# Where y = [x; s; t]
# --------------------------
# lambdaInit = zeros(size(fObj(x0)))
# rhoInit = 1

@doc raw"""
    augLag(obj::objectiveFunc, cM::constraintManager, rho::Float64)

This is a mutable struct that stores an SOCP Augmented Lagrangian.

For a single SOCP constraint, the lagrangian is
```math
φ(y) = f(y) + (ρ/2) ||c(y)||_2^2 + λ c(y)
```
```math
φ(y) = f(y) + (ρ/2) c(y)'c(y)    + λ c(y)
```

Where `y` is represented a primal vector struct.

Generally, we initialize `ρ` to 1.

See also: [`evalAL`](@ref), [`evalGradAL`](@ref), [`evalHessAl`](@ref)
"""
mutable struct augLag
    obj::objectiveFunc
    cM::constraintManager
    rho::Float64
end


function getPrimals(al::augLag, primalDualVec)
    if typeof(al.cM) == constraintManager_Dynamics
        # println("Adjusting Primal Vec")
        affLambdaSize = size(al.cM.affineLambdaList, 1)
        primalSize = size(primalDualVec, 1) - affLambdaSize
        primals = primalDualVec[1:primalSize]
    else
        # There are no duals
        primals = primalDualVec
    end

    return primals
end


"""
    evalAL(alQP::augLag, y)

Evaluates the Augmented Lagrangian with the primals `y`. If linear dynamics are
present, it will special case the dynamics.

returns a real number
"""
function evalAL(alQP::augLag, y)
    # φ(y) = f(x) + (ρ/2) c(y)'c(y)    + λ c(y)
    # φ(y) = [(1/2) xT Q x + cT x] + (ρ/2) (c(y))'(c(y)) + λ (c(y))

    primals = getPrimals(alQP, y)

    fCurr = fObjQP(alQP.obj, primals)[1]
    cCurr = evalConstraints(alQP.cM, primals, alQP.rho)

    # println("f(x) = $fCurr")
    # println(" with size -> $(size(fCurr))")
    # println("c(x) = $cCurr")
    # println(" with size -> $(size(cCurr))")

    if typeof(alQP.cM) == constraintManager_Dynamics
        # φ(y) + λ d(y)
        duals = alQP.cM.affineLambdaList
        cAff = evalAffineEq(alQP.cM, primals)
        # println("Duals (λ) = $(size(duals))")
        # println("Ax - b = $(size(cAff))")
        cAffVal = duals'cAff
        # println("λ'(Ax-b) = $(size(cAffVal)) = $cAffVal")
        # println("-> $(cAffVal[1])")
        cCurr += cAffVal[1]
    end

    return fCurr + cCurr
end

"""
    evalGradAL(alQP::augLag, y)

Evaluates the gradient of the augmented lagrangian with the primals `y`. The
format is adjusted depending on the type of constraint manager.

For a base constraint manager (all constraints in the augmented lagrangian),
the gradient is just ∇φ = ∇f + Σ ∇ci.

For a constraint manager that separates the equality constraints, the gradient
also includes the derivatives with respect to the duals of the equality
constraints. So ∇φ = [∇f + Σ ∇ci; Ax - b]

returns a vector of size `y`
"""
function evalGradAL(alQP::augLag, y)
    # ∇φ(y) = ∇f(x) + J(c(y))'(ρ c(y) + λ)
    # y = [x; s; t]
    # ∇φ(y) is (n+m+1)x1

    primals = getPrimals(alQP, y)

    gradfCurr = dfdxQP(alQP.obj, primals)
    gradCCurr = evalGradConstraints(alQP.cM, primals, alQP.rho)

    gradPhiPrimal = gradfCurr + gradCCurr

    if typeof(alQP.cM) == constraintManager_Base
        return gradPhiPrimal
    elseif typeof(alQP.cM) == constraintManager_Dynamics
        duals = alQP.cM.affineLambdaList
        Aaff = affineBlock(alQP.cM)
        return [gradPhiPrimal + Aaff'duals; evalAffineEq(alQP.cM, primals)]
    end
end

"""
    evalHessAl(alQP::augLag, y)

Evaluates the Hessian of the Augmented Lagrangian of an SOCP with the
primals `y`. The format is adjusted depending on the type of constraint manager.

For a base constraint manager (all constraints in the augmented lagrangian),
the hessian is just Hφ = Hf + Σ Hci.

For a constraint manager that separates the equality constraints, the gradient
also includes the derivatives with respect to the duals of the equality
constraints. So Hφ = [(Hf + Σ Hci)     A'; A    0]

returns a Symmetric Matrix of size `y`×`y` or `y + b`×`y + b`
"""
function evalHessAl(alQP::augLag, y)

    primals = getPrimals(alQP, y)

    hessf = hessQP(alQP.obj)
    hessC = evalHessConstraints(alQP.cM, primals, alQP.rho)

    hessPhiPrimals = hessf + hessC

    if typeof(alQP.cM) == constraintManager_Base
        return hessPhiPrimals
    elseif typeof(alQP.cM) == constraintManager_Dynamics
        Aaff = affineBlock(alQP.cM)
        blockSize = size(Aaff, 1)
        return [hessPhiPrimals Aaff'; Aaff spzeros(blockSize, blockSize)]
    end
end

"""
    calcNormGradResiduals(alQP::augLag, yList)

Calculates the 2-norm of the residuals of an augmented lagrangian at each of
the inputs `y`. The 2-norm in safe in that it is bounded below (e.g. by
10^(-20)) and can be plotted on a semilog or log-log plot.

See also: [`safeNorm`](@ref)
"""
function calcNormGradResiduals(alQP::augLag, yList)
    #=
    Calculate the residuals where the AL is the merit function.
    Returns the norm of the gradient of the AL at each point in xArr
    =#

    primals = getPrimals(alQP, yList[1])
    pSize = size(primals, 1)
    resArr = [evalGradAL(alQP, y[1:pSize]) for y in yList]

    return safeNorm(resArr)
end

"""
    calcALArr(alQP::augLag, yList::Array{SOCP_primals, 1})

Evaluates the augmented lagrangian at each input `y`.
"""
function calcALArr(alQP::augLag, yList)
    #=
    Calculate the value of the Augmented Lagrangian at each point
    =#
    return [evalAL(alQP, y) for y in yList]
end
