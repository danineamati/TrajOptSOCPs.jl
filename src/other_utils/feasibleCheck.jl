

# This simply holds the feasibility checks
function isFeasiblePolyHedron(A, b, x)
    #=
    Checks if all polyhedron checks are satisfied
    Namely: Ax ≤ b
    =#

    state = A * x - b

    if maximum(state) > 0
        return false
    end

    return true
end
