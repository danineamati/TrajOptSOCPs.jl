# Implements trust region code

using LinearAlgebra, SparseArrays


@doc raw"""
    dampingInitialization(B, g, delta, epsilon, a = 0.5)

Finds an appropriate damping parameter to start the search.

Specifically,
```math
λ ∈ [0, ||B|| + (1 + ϵ) \frac{||g||}{Δ}]
```

## Arguments
- `B`: The (approximate or true) Hessian Matrix
- `g`: The (approximate or true) gradient vector
- `epsilon`: Parameter (> 0) weighes the gradient versus the hessian in
determining the damping ratio size
- `a`: Parameter (> 0) How quickly to damp the Hessian

See also: [`findDamping`](@ref)
"""
function dampingInitialization(B, g, delta, epsilon, IPD, a = 0.5,
                                    basedamping = 1e-6)

    dampingMax = (norm(B) + ((1 + epsilon) * norm(g) / delta))

    if isposdef(B)
        damping = basedamping
        Bdamped = B + damping * IPD
        cholBdamp = cholesky(Bdamped, check = false)

        if issuccess(cholBdamp)
            return damping, dampingMax, cholBdamp
        end
    end

    while true
        damping = a * dampingMax

        Bdamped = B + damping * IPD

        # Get the Cholesky decomposition without throwing an exception
        luB = lu(Bdamped, check = false)

        # Checks positive definite and non-singular
        if issuccess(luB)
            return damping, dampingMax, Bdamped
        else
            # If the it fails, increase damping factor. Converge to a -> 1
            a = (a + 1)/ 2
            println("More Damping...")

            # Once "a" gets close to 1, might as well stop
            if a > 0.99
                return dampingMax, dampingMax, Bdamped
            end
        end
    end
end

"""
    findDamping(B, g, delta = 1, gamma = 1.5, epsilon = 0.5,
                      condNumMax = 1e7, a = 1e-6, verbose = false)

Implements a search for the appropriate trust region.

## Arguments:
- `B`: The (approximate or true) Hessian Matrix
- `g`: The (approximate or true) gradient vector
- `delta`: is the initial trust region size
- `condNumMax`: The max condition number on the damped Hessian

for `gamma`, `epsilon`, and `a`, see [`dampingInitialization`](@ref)

returns the damping factor `damping`, the newton step `dk`, and the damped
hessian represented as a Cholesky Decomposition `rCho`

This corresponds to algorithm 2.6 in Nocedal et Yuan (1998)

"""
function findDamping(B, g, dualSize = 0, delta = 1, gamma = 1.5, epsilon = 0.5,
                        condNumMax = 1e7, a = 1e-6, verbose = false)

    # For Primal-Dual the effective Identity is [I 0; 0 -I] (albeit with
    # correct dimensions)
    getIPD(p, d) = Diagonal([ones(p); -ones(d)])
    IPD = getIPD(size(g, 1) - dualSize, dualSize)

    # Initialize the damping factor and find the initial damped hessian
    # represented as a Cholesky Decomposition (rCho)
    # damping, dampingMax, rCho = dampingInitialization(B, g, delta, epsilon,
    #                                                   IPD, a)
    damping, dampingMax, BDamp = dampingInitialization(B, g, delta, epsilon,
                                                      IPD, a)

    # Use the Cholesky decomposition to more easily solve the problem.
    # LCho = sparse(rCho.L)
    # dk = - LCho' \ (LCho \ g)
    dk = - BDamp \ g

    if verbose
        println("Damping Start = $damping and dampingMax = $dampingMax")
        println("dk = $dk")
    end

    counter = 1

    # Have to use Array to convert briefly to a dense matrix in order
    # to determine the condition number.
    # while (norm(dk) > delta) && (cond(Array(LCho * LCho')) > condNumMax)
    while (norm(dk) > delta) && (cond(Array(BDamp)) > condNumMax)
        if verbose
            println("Increasing Damping")
            println("Cholesky Factorization:")
            display(rCho)
            println("dk")
            display(dk)
        end

        LCho = lu(BDamp).L

        qk = LCho \ dk # Problematic

        updateNumerator = norm(dk)^2 * (gamma * norm(dk) - delta)
        updateDenominator = norm(qk)^2 * delta

        damping = min(damping + updateNumerator / updateDenominator, dampingMax)

        if verbose
            println("qk = $qk")
            println("New Damping: $damping")
        end

        # rCho = lu(B + damping * IPD, check = false)
        BDamp = B + damping * IPD
        luB = lu(BDamp, check = false)

        if issuccess(luB)
            dk = - BDamp \ g
            counter += 1

            if damping == dampingMax
                break
            end
        else
            # Need to really increase the damping factor
            damping *= 10
        end


        if counter >= 10
            break
        end
    end

    println("Condition Number: $(cond(Array(BDamp)))")

    return damping, dk, BDamp

end
