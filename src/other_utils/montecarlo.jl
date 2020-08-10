

# This function takes care of the Monte Carlo aspects of the tests
using Distributions

function getPoint(minVec, maxVec)
    #=
    Returns a random point in the rectangle specified by the min and max
    vectors

    In this case, a uniform distribution is used. In future work, a different
    distribution could be used without loss of generality.
    =#
    pt = zeros(size(minVec))
    for dimen in 1:size(minVec, 1)
        pt[dimen] = rand(Uniform(minVec[dimen], maxVec[dimen]))
    end
    return pt
end

function montecarlo(minVec, maxVec, numPoints = 10, funcCheck = x -> true)
    #=
    This function takes a given range and randomly selects a point in that
    range. It then checks if that point passes a check before adding it.
    The result is a list of random points.

    In general, the check function is for the feasible set.
    Ex: x -> isFeasiblePolyHedron(A1, b1, x)

    If the check function is left as default, all points pass and are added
    to the list.
    =#

    listPoints = []

    while size(listPoints, 1) < numPoints
        newpt = getPoint(minVec, maxVec)

        if funcCheck(newpt)
            push!(listPoints, newpt)
        end
    end

    return listPoints
end

function inBox(x2D)
    if (1 > x2D[1] > 0) && (1 > x2D[2] > 0)
        return true
    end
    return false
end
