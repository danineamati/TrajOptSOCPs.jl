# Linear Quadratic Regulator Objective Functions

using SparseArrays

include("QP_Linear_objectives.jl")

abstract type objectiveLQR_abstract <: objectiveFunc end

"""
    LQR_simple

A simple LQR struct where `Q` and `R` are constant for all `x` and `u`.
"""
struct LQR_simple <: objectiveLQR_abstract
    Q
    R
end

"""
    LQR_QP

An LQR what acts like a QP objective function.
"""
struct LQR_QP <: objectiveLQR_abstract
    QR_Full
end

"""
    LQR_QP

An LQR what acts like a QP objective function, but with a reference
trajectory.
"""
struct LQR_QP_Referenced <: objectiveLQR_abstract
    lqr_qp::LQR_QP
    xRef
end


#####################################
#  LQR -> ObjectiveQP Constructors  #
#####################################
@doc raw"""
    makeLQR_TrajSimple(Q, R, NSteps::Int64)

Pseudo-constructor for an LQR objective function (which is really an
[`objectiveQ`](@ref) objective function).

An LQR struct for a trajectory is one that holds the sparse matrix matching:

```math
B_{QR} = [Q 0 0 0; 0 R 0 0; 0 0 Q 0; 0 0 0 R]
```

Where size corresponds to the size of the trajectory (or trajectory horizon).

```math
f(y) = \frac{1}{2} \ y^{\top} * B_{QR} * y
```

See also: [`LQR_simple`](@ref) and [`LQR_QP`](@ref)
"""
function makeLQR_TrajSimple(Q, R, NSteps::Int64)
    QSize = size(Q, 1)
    RSize = size(R, 1)
    totSize = (QSize + RSize) * NSteps + QSize

    # Initialize and empty sparse array
    QRFull = spzeros(totSize, totSize)

    rStart = 1

    for k in 1:NSteps
        rEnd = rStart + QSize - 1
        QRFull[rStart:rEnd, rStart:rEnd] = Q

        rStart += QSize
        rEnd = rStart + RSize - 1

        QRFull[rStart:rEnd, rStart:rEnd] = R

        rStart += RSize
    end
    rEnd = rStart + QSize - 1
    QRFull[rStart:rEnd, rStart:rEnd] = Q

    return LQR_QP(QRFull)
end

function makeLQR_TrajSimple(lqr::LQR_simple, NSteps::Int64)
    return makeLQR_TrajSimple(lqr.Q, lqr.R, NSteps)
end

function makeLQR_TrajReferenced(Q, R, NSteps::Int64, xRef)
    baseLQR = makeLQR_TrajSimple(Q, R, NSteps)
    return LQR_QP_Referenced(baseLQR, xRef)
end

function makeLQR_TrajReferenced(lqr::LQR_simple, NSteps::Int64, xRef)
    baseLQR = makeLQR_TrajSimple(lqr.Q, lqr.R, NSteps)
    return LQR_QP_Referenced(baseLQR, xRef)
end


#######################################
###       Evaluate An LQR QP        ###
#######################################

"""
    fObjQP(lqr::LQR_QP, x)

Evaluates a full QP-like LQR. See [`LQR_QP`](@ref).
"""
function fObjQP(lqr::LQR_QP, x)
    return (1/2) * x' * (lqr.QR_Full) * x
end

"""
    fObjQP(lqr::LQR_QP, x, xRef)

Evaluates a full QP-like LQR relative to a reference. (i.e. x -> (x - xRef))

See [`LQR_QP`](@ref).
"""
function fObjQP(lqr::LQR_QP, x, xRef)
    xDiff = x - xRef
    return fObjQP(lqr, xDiff)
end

"""
    fObjQP(lqrRef::LQR_QP_Referenced, x)

Evaluates a full QP-like LQR relative to a reference. (i.e. x -> (x - xRef))

See [`LQR_QP`](@ref), [`LQR_QP_Referenced`](@ref)
"""
function fObjQP(lqrRef::LQR_QP_Referenced, x)
    # println("x = $(size(x))")
    # println("xRef = $(size(lqrRef.xRef))")
    xDiff = x - lqrRef.xRef
    return fObjQP(lqrRef.lqr_qp, xDiff)
end

#######################################
###      Gradients of LQR QP        ###
#######################################
"""
    dfdxQP(lqr::LQR_QP, x)

Evaluates the derivative of a full QP-like LQR at input `x`.

See [`LQR_QP`](@ref), [`dfdxQP`](@ref)
"""
function dfdxQP(lqr::LQR_QP, x)
    return (lqr.QR_Full) * x
end

"""
    dfdxQP(lqr::LQR_QP, x, xRef)

Evaluates the derivative of a full QP-like LQR at input `x` with respect to
a reference `xRef`.

See [`LQR_QP`](@ref), [`dfdxQP`](@ref)
"""
function dfdxQP(lqr::LQR_QP, x, xRef)
    xDiff = x - xRef
    return dfdxQP(lqr, xDiff)
end

"""
    dfdxQP(lqrRef::LQR_QP_Referenced, x)

Evaluates the derivative of a full QP-like LQR at input `x` with respect to
a reference `xRef`.

See [`LQR_QP_Referenced`](@ref), [`LQR_QP`](@ref), [`dfdxQP`](@ref)
"""
function dfdxQP(lqrRef::LQR_QP_Referenced, x)
    return dfdxQP(lqrRef.lqr_qp, x, lqrRef.xRef)
end

#######################################
###      Hessians of LQR QP         ###
#######################################
"""
    hessQP(lqr::LQR_QP)

Evaluates the hessian of a full QP-like LQR, which is simple the base matrix.

See [`LQR_QP`](@ref)
"""
function hessQP(lqr::LQR_QP)
    return lqr.QR_Full
end

"""
    hessQP(lqr::LQR_QP)

Evaluates the hessian of a full QP-like LQR, which is simple the base matrix.

See [`LQR_QP`](@ref)
"""
function hessQP(lqrRef::LQR_QP_Referenced)
    return hessQP(lqrRef.lqr_qp)
end
