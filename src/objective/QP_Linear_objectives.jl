# Holds the possible QP-like objective structs

#############################################################
###  Structs for Quadratic and Linear Objective Functions ###
#############################################################

abstract type objectiveFunc end
abstract type objectiveQP_abstract <: objectiveFunc end

using SparseArrays



@doc raw"""
    objectiveQP(Q, p)

Quadratic (and Linear) Objective/Cost Function
----

```math
\frac{1}{2} x^{⊤} Q x + p^{⊤} x
```

`Q` is an nxn Symmetric Matrix (that is positive definite for convex)

`p` is an nx1 vector

See also: [`objectiveQP`](@ref), [`fObjQP`](@ref), [`dfdxQP`](@ref)
"""
struct objectiveQP <: objectiveQP_abstract
    "Positive Definite Symmetric Matrix"
    Q
    p
end

@doc raw"""
    objectiveQ(Q)

Quadratic (only) Objective/Cost Function
----

```math
\frac{1}{2} x^{⊤} Q x
```

`Q` is an nxn Symmetric Matrix (that is positive definite for convex)

See also: [`objectiveQP`](@ref) [`fObjQP`](@ref), [`dfdxQP`](@ref)
"""
struct objectiveQ <: objectiveQP_abstract
    "Positive Definite Symmetric Matrix"
    Q
end

@doc raw"""
    objectiveP(P)

Linear (only) Objective/Cost Function
----

```math
p^{⊤} x
```

`p` is an nx1 vector

See also: [`objectiveQP`](@ref), [`fObjQ`](@ref), [`dfdxQ`](@ref)
"""
struct objectiveP <: objectiveQP_abstract
    "Positive Definite Symmetric Matrix"
    p
end

#############################################################
###       Evaluate Quadratic and Linear Objectives        ###
#############################################################


"""
    fObjQP(qp::objectiveQP, x)

Evaluates a quadratic function at the input `x`. The input `x` should be a
*column* vector (i.e. `x = [4; 3]`)
"""
function fObjQP(qp::objectiveQP, x)
    return (1/2) * x' * (qp.Q) * x + (qp.p)' * x
end

"""
    fObjQP(q::objectiveQ, x)

Evaluates a quadratic (only) function at the input `x`. The input `x` should
be a *column* vector (i.e. `x = [4; 3]`)
"""
function fObjQP(q::objectiveQ, x)
    return (1/2) * x' * (q.Q) * x
end

"""
    fObjQP(p::objectiveP, x)

Evaluates a linear function at the input `x`. The input `x` should be a
*column* vector (i.e. `x = [4; 3]`)
"""
function fObjQP(p::objectiveP, x)
    return (p.p)' * x
end

#############################################################
###     Gradients for Quadratic and Linear Objectives     ###
#############################################################


@doc raw"""
    dfdxQP(qp::objectiveQP, x)

Evaluates the derivative of a quadratic function at input `x`. The input `x`
should be a *column* vector (i.e. `x = [4; 3]`)

```math
\frac{d}{dx} (\frac{1}{2} x^{⊤} Q x + p^{⊤} x) = Qx + p
```
"""
function dfdxQP(qp::objectiveQP, x)
    return (qp.Q) * x + qp.p
end

@doc raw"""
    dfdxQ(q::objectiveQ, x)

Evaluates the derivative of a quadratic (only) function at input `x`. The
input `x` should be a *column* vector (i.e. `x = [4; 3]`)

```math
\frac{d}{dx} (\frac{1}{2} x^{⊤} Q x) = Qx
```
"""
function dfdxQP(q::objectiveQ, x)
    return (q.Q) * x
end

@doc raw"""
    dfdxQP(p::objectiveP, x = 0)

Evaluates the derivative of a linear function. The `x` is unused, but is kept
to maintain the same function header.

```math
\frac{d}{dx} (p^{⊤} x) = p
```
"""
function dfdxQP(p::objectiveP, x = 0)
    return p.p
end


#############################################################
###      Hessians for Quadratic and Linear Objectives     ###
#############################################################


@doc raw"""
    hessQP(qp::objectiveQP)

Evaluates the hessian of a quadratic function.

```math
H(\frac{1}{2} x^{⊤} Q x + p^{⊤} x) = Q
```

Note that Q is already square, symmetric, and (generally) positive definite.
"""
function hessQP(qp::objectiveQP)
    return qp.Q
end

@doc raw"""
    hessQP(q::objectiveQ)

Evaluates the hessian of a quadratic (only) function.

```math
H(\frac{1}{2} x^{⊤} Q x) = Q
```

Note that Q is already square, symmetric, and (generally) positive definite.
"""
function hessQP(q::objectiveQ)
    return q.Q
end

@doc raw"""
    hessQP(q::objectiveQ)

Evaluates the hessian of a linear (only) function. (which is zero)

```math
H(p^{⊤} x) = 0
```
"""
function hessQP(p::objectiveP)
    n = size(p.p, 1)
    return spzeros(n, n)
end

#############################################################
