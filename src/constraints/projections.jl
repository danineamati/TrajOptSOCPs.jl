
using LinearAlgebra

# Covers how to consider a square violation
"""
    sqrDistEuclid(pt, projFunc, p = 2)

Provide a point `pt` and a function `pt -> projFunc(pt)` and this function will
return the squared distance to that point. A `p = 2` norm is assumed, but that
can be adjusted if needed.

Primarily used to visualize square penalty constraints.
"""
function sqrDistEuclid(pt, projFunc, p = 2)
    #=
    Choose a given projection method and determine the squared distance
    at this point.
    =#
    return norm(pt - projFunc(pt), p)^2
end



# Covers Projections

# First, Lines from vectors at the origin
@doc raw"""
    projLine(pt, vec)

Project a point `pt` onto a line described by `vec`. `pt` should have the same
dimensionality as `vec`.

```math
Proj_v(p) = v \ \frac{p ⋅ v}{v ⋅ v}
```

See also: [`orthoProjLine`](@ref)
"""
function projLine(pt, vec)
    #=
    This function projects a point pt onto a line described by vec

    proj_v(p) = v * (p'v)/(v'v)
    =#
    return vec * (pt'vec) / (vec'vec)
end

@doc raw"""
    orthoProjLine(pt, vec)

Project a point `pt` onto the orthogonal of a line described  by `vec`. `pt`
should have the same dimensionality as `vec`.

```math
Proj_{v_⟂}(p) = p - Proj_v(p)
```

See also: [`projLine`](@ref)
"""
function orthoProjLine(pt, vec)
    #=
    This function projects a point pt onto the orthogonal of a line described
    by vec

    proj_v(p) = v * (p'v)/(v'v)
    proj_v⟂(p) = pt - v * (p'v)/(v'v)
    =#
    return pt - projLine(pt, vec)
end

# Second, Positive Orthant
@doc raw"""
    projPosOrth(pt)

Project a point `pt` onto the positive orthant. This is most commonly used for
constraints of the form `τ ≥ 0`.

```math
Proj_+(p) = \max.(p, 0)
```

"""
function projPosOrth(pt)
    #=
    This function is for a constraint of the form "τ ≥ 0."
    When τ_i ≥ 0, the variable τ is in the positive orthant, thus, the
    projection is simply τ_i.
    When τ_i < 0, the variable τ is not in the positive orthant, thus, the
    projection is distance to the positive orthant
    =#

    return max.(pt, 0)
end

# Third, Affine Case
@doc raw"""
    projAffineEq(a, b, x)

Project a point `x` onto an n-dimensional hyperplane. Equivalently, this is a
single affine constraint of the form `a'x = b`. Note that `a` and `x` should
have the same dimensionality and `b` is a real number. This function does
NOT handle matrix versions of the form `Ax = b`. To handle cases of this form,
split `A` and `b` by row and pass the same `x`. In particular:

Turn `Ax = b` into `A[1, :]'x = b1`, `A[2, :]'x = b2`, etc.

Important, there is a singularity at `a = 0` corresponding to a hyperplane
independent of `x`! Please avoid this.

```math
Proj_{Ax = b}(x) = x + \frac{ba}{a ⋅ a} - a \ \frac{a ⋅ x}{a ⋅ a}
```

See also [`projAffineIneq`](@ref)
"""
function projAffineEq(a, b, x)
    #=
    Projection for Affine Equality Case (a'x = b)

    We can show that the projection is given by
    proj(x) = x + (b - a'x) * a / a'a
            = x + a * [b] / a'a - a * (a'x) / a'a
    =#
    if a'a == 0
        println("Divide by zero error in projection, returning zero.")
        return 0
    end
    t1 = a * b / a'a
    # println("T1 = $t1")
    t2 = a * (a'x) / a'a
    # println("T2 = $t2")
    return x + t1 - t2
end

"""
    projAffineIneq(a, b, x)

Project a point `x` onto an n-dimensional half space. Equivalently, this is a
single affine constraint of the form `a'x ≤ b`. This function is entirely based
on [`projAffineEq`](@ref).

Please see [`projAffineEq`](@ref) for details. When the constraint is satisfied,
the projection is simply `x` itself.
"""
function projAffineIneq(a, b, x)
    #=
    Projection for Affine Inequality Case (a'x ≤ b)

    We can show that the projection is given by
    projAffineEq(x) for a'x > b
    x               for a'x ≤ b
    =#
    # println("a = $(size(a))")
    # display(a)
    # println("x = $(size(x))")
    if (a'x)[1] > b
        return projAffineEq(a, b, x)
    end
    return x
end

@doc raw"""
    projSecondOrderCone(v, s, filled = true, p = 2)

Project `(v, s)` onto the second order cone described by
```math
C = \{(x, t) ∈ R^{n + 1} \ | \ ||x||_2 ≤ t\}
```
if the cone is filled (i.e. `filled = true`). Otherwise, project onto the
boundary of the second order cone described by
```math
C = \{(x, t) ∈ R^{n + 1} \ | \ ||x||_2 = t\}
```

The parameter `p` would allow you to select a different cone, however, for now
only the second order cone (`p = 2`) is functional.

---
Three projection cases:
1. The point is below the tip. This returns 0.
2. The point is interior to the filled cone. This returns `(v, s)` in the
filled case and projects onto the boundary in the unfilled case.
3. The point is exterior to the cone, but not below the tip. This projects
onto the boundary of the cone.

```math
Proj_{Cone \ Border}(v, s) = \frac{1}{2} (1 + \frac{s}{||v||}) \ [v, ||v||]
```

"""
function projSecondOrderCone(v, s, filled = true, p = 2)
    #=
    Projection for Second-Order Cone (AKA quadratic cone or the Lorentz cone)

    The second-order cone is C = {(x, t) ∈ R^n+1 | ||x||2 ≤ t}. Using the
    2-norm. Projection onto it is given by

    proj(v, s) =
    0                                   for ||v|| ≤ -s  (Below the tip)
    (v, s)                              for ||v|| ≤ s   (In the cone)
    (1/2)(1 + s/||v||)(v, ||v||)        for ||v|| ≥ |s| (Onto the cone)

    Note that (|s| = absolute value of s)
    =#

    if norm(v, p) ≤ -s
        # println("Below Tip")
        return zeros(size([v; s]))
    elseif norm(v, p) ≤ s
        # println("Inside already")
        if filled
            return [v; s]
        end
        if !filled
            return (1/2) * (1 + s / norm(v, p)) * [v; norm(v, p)]
        end
    elseif norm(v, p) ≥ abs(s)
        # println("Outside")
        return (1/2) * (1 + s / norm(v, p)) * [v; norm(v, p)]
    end

    # The above cases should be all inclusive, but as a safety I put an error.
    error("Second Order Cone Conditions ERROR")

end
