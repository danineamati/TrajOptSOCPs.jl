# KKT QP with Equality Constraints Solver

include("..\\auglag\\auglag-core.jl")

using SparseArrays


@doc raw"""
    solveQP_AffineEq

We have a problem of the form:

```math
\underset{x}{\text{minimize}} \quad \frac{1}{2} x^⊤Q x + p^⊤x
```
```math
\text{subject to} \quad Ax = b
```

This yields the KKT system:

```math
[Q A^⊤; A 0] \ [x; λ] = [-p; b]
```

This can be solved with LU Factorization, but we will let Julia's "\" function
handle that.

"""
function solveQP_AffineEq(Q, p, A, b)
    hMat = [Q A'; A spzeros(size(A, 1), size(A, 1))]
    dMat = [-p; b]

    return hMat \ Array(dMat)
end

function solveQP_AffineEq(obj::objectiveQP, constraints::AL_AffineEquality)
    Q = obj.Q
    p = obj.p

    A = constraints.A
    b = constraints.b

    return solveQP_AffineEq(Q, p, A, b)
end
