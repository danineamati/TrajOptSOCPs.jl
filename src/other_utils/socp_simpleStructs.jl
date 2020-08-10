# SOCP Primal Structs for Simple SOCP problems

include("..\\constraints\\constraints.jl")



"""
    SOCP_primals(x, s, t)

Struct to store the primal variables of the form [x; s; t]

`x` and `s` are mx1 and nx1 vectors.
`t` is a real number

See also: [`primalVec`](@ref) and [`primalStruct`](@ref)
"""
mutable struct SOCP_primals
    x
    s
    t
end

"""
    primalVec(y::SOCP_primals)

Converts [`SOCP_primals`](@ref) struct to [x; s; t] vector.
"""
function primalVec(y::SOCP_primals)
    return [y.x; y.s; y.t]
end

"""
    primalStruct(v, xSize::Int64, sSize::Int64, tSize::Int64)

Converts [x; s; t] vector to [`SOCP_primals`](@ref) struct
"""
function primalStruct(v, xSize::Int64, sSize::Int64, tSize::Int64)
    return SOCP_primals(v[1:xSize], v[xSize+1:xSize+sSize], v[end])
end


"""
    getXVals(yList::Array{SOCP_primals, 1})

Returns the "x" values of and an array of [`SOCP_primals`](@ref).
Used primarily with plotting.
"""
function getXVals(yList::Array{SOCP_primals, 1})
    xList = []
    for xst in yList
        push!(xList, xst.x)
    end
    return xList
end


@doc raw"""
    getViolation(yList::Array{SOCP_primals, 1}, c::AL_coneSlack)

returns how much each SOCP constraint is violated
(in 3 mulitdimensional arrays).

Recall that an SOCP is generally
```math
||Ax-b|| ≤ c^⊤ x - d
```

Using slack variables, we have

```math
||s|| - t ≤ 0
```

```math
(Ax - b) - s = 0
```

```math
(c^⊤ x - d) - t = 0
```

This function returns the constraint violation for each of the above 3 in 3
separate arrays.
"""
function getViolation(yList::Array{SOCP_primals, 1}, c::AL_coneSlack)
    coneViolation = []
    affineViolation = []
    lastViolation = []
    for y in yList
        vio = getNormToProjVals(c, y.x, y.s, y.t)
        push!(coneViolation, vio[1])
        push!(affineViolation, vio[2:1+size(y.s, 1)])
        push!(lastViolation, vio[end])
    end
    return coneViolation, affineViolation, lastViolation
end
