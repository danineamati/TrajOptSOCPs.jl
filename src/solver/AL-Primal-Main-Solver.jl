#=

Augmented Lagrangian solver.

It is not abstracted as to handle multiple constraints

=#


using LinearAlgebra

include("..\\auglag\\auglag-core.jl")
include("..\\auglag\\solverParams.jl")

# include("..\\constraints\\constraintManager.jl")

include("backtrackLineSearch.jl")
include("trustRegion.jl")


@doc raw"""
    newtonStepALP(y0::SOCP_primals, al::augLagQP_2Cone, delta = 1,
                      gamma = 1.5, epsilon = 0.5, condNumMax = 1e5, a = 1e-6)

Evaluates the Newton Step for an Augmented Lagrangian of an SOCP.

```math
y ← y - [H(φ(y)) + λI]^{-1} ∇φ(y)
```

This function internally calculates the step `dk` using Cholesky Decomposition.
```math
dk = [H(φ(y)) + λI] \ ∇φ(y)
```

## Core Arguments
- `y0::SOCP_primals`: The primal vectors of the SOCP.
- `al::augLagQP_2Cone`: The SOCP Augmented Lagrangian.
See: [`SOCP_primals`](@ref), [`augLagQP_2Cone`](@ref)

## Trust Region and Regularization Arguements
- `delta`: Current trust region size
- `gamma`: Parameter (> 1) for how aggressive to increase the trust region size
- `epsilon`: Parameter (> 0) weighes the gradient versus the hessian in
determining the damping ratio size
- `condNumMax`: The max condition number on the damped Hessian
- `a`: Parameter (> 0) How quickly to damp the Hessian

See: [`findDamping`](@ref) for more information


## Return Values
returns `dk`, `phiD`, `rCho`, `damping` which are
- `dk`: the Newton step
- `phiD`: the residual (``∇φ(y)``)
- `rCho`: the Cholesky decomposition of `H(φ(y)) + λI`
- `damping`: the damping parameter (`λ` above)

See: [`findDamping`](@ref) for more information

"""
function newtonStepALP(y0, al::augLag, delta = 1,
                    gamma = 1.5, epsilon = 0.5, condNumMax = 1e5, a = 1e-6)
    hess = evalHessAl(al, y0)
    phiD = evalGradAL(al, y0)

    if typeof(al.cM) == constraintManager_Dynamics
        dualSize = size(al.cM.affineLambdaList, 1)
    else
        dualSize = 0
    end

    damping, dk, rCho = findDamping(hess, phiD, dualSize, delta, gamma, epsilon,
                                            condNumMax, a)

    if norm(dk) > delta
        dk = (delta / norm(dk)) * dk
    end

    return dk, phiD, rCho, damping

end

"""
    newtonTRLS_ALP(y0::SOCP_primals, al::augLagQP_2Cone,
                            sp::solverParams, verbose = false)

Performs a maximum on `N` newton steps as specified by [`solverParams`](@ref).
May return early if the residual is lower than the tolerance specified by
[`solverParams`](@ref).

The function operates on an augmented lagrangian (AL) with primals (P) only
(as opposed to a Primal-Dual set-up) for second-order cone programs (SOCP).
The function uses a combined trust region (TR) and line search (LS) approach
as specified in Nocedal et Yuan (1998).

returns (yNewtStates, residNewt) which are:
- `yNewtStates`: An array with one entry of the primal vector per Newton step
- `residNewt`: An array with one entry of the residual vector per Newton step

See also: [`augLagQP_2Cone`](@ref), [`SOCP_primals`](@ref)
"""
function newtonTRLS_ALP(y0, al::augLag, sp::solverParams, verbose = false)
    yNewtStates = []
    residNewt = []
    # push!(yNewtStates, y0)
    yCurr = y0

    trustDelta = sp.trSizeStart

    # For printing
    early = false

    lineSearchObj(v) = evalAL(al, v)
    lineSearchdfdx(v) = evalGradAL(al, v)

    # Now, we run through the iterations
    for i in 1:(sp.maxNewtonSteps)

        if verbose
            println("Currently at $yCurr")
        end

        currentObjVal = lineSearchObj(yCurr)

        if verbose
            println("ϕ(y) = $currentObjVal")
            println("∇ϕ(y) = $(lineSearchdfdx(yCurr))")
            cCurr = evalConstraints(al.cM, yCurr, al.rho)
            println("Constraints: $cCurr")
        end

        # Take a Newton Step
        # Negative sign addressed above
        # (dk, residual, rCho, damp) = newtonStepALP(yCurr, al, trustDelta)
        (dk, residual, BDamp, damp) = newtonStepALP(yCurr, al, trustDelta)
        push!(residNewt, residual)
        # LCho = sparse(rCho.L)

        if true
            # println("\nNewton Direction: $dk")
            println("Damping of ($damp) and trust size of ($trustDelta)")
            # println("AL: $al")
        end

        # Determine if the trust region is sufficient
        y0New = yCurr + dk
        baseObjVal = lineSearchObj(y0New)
        println("Objective from $currentObjVal → $baseObjVal")
        if baseObjVal ≤ currentObjVal
            # Trust region was a success
            # Step 4 of algorithm 3.1 in the Nocedal et Yuan paper

            if true
                println("Trust Region Success")
            end

            # fObjApprox(d) = residual'd + (1/2) * d'*(LCho * LCho')*d
            fObjApprox(d) = residual'd + (1/2) * d'*(BDamp)*d
            println("fObjApprox(d) = $(fObjApprox(dk))")
            rho = (currentObjVal - baseObjVal) / (-fObjApprox(dk)[1])
            roundingErrorTol = UInt8(round(-log10(sp.rTol)) / 2)

            println("Rho = $rho")
            roundedDkNorm = round(norm(dk), digits = roundingErrorTol)
            roundedTrustDelta = round(trustDelta, digits = roundingErrorTol)

            print("||dk|| = $(norm(dk)) → $roundedDkNorm")
            print(" vs Δ = $trustDelta → $roundedTrustDelta")
            println(" at $roundingErrorTol digits")

            if rho ≥ sp.trc2 && roundedDkNorm < roundedTrustDelta
                # ρ ≥ c2 AND ||dk|| < Δ
                # Condition 1: No change
                # trustDelta = trustDelta # Probably not needed
                println("Trust Region Size unchanged. Still at $trustDelta")
            elseif rho < sp.trc2
                # ρ < c2, but ||dk|| ≤ Δ
                trustDelta = (sp.trc3 * norm(dk) +
                                            sp.trc4 * trustDelta) / 2
                println("Trust Region Size decreased. Now $trustDelta")
            else
                # Implies that ρ ≥ c2 AND ||dk|| = Δ
                trustDelta = sp.trc1 * trustDelta
                println("Trust Region Size increased. Now $trustDelta")
            end


        else
            # Trust region failed.

            if true
                print("Trust Region Failed - ")
                println("Trying Line Search.")
            end

            # println("Attempting Line Search: ")
            lc = lineSearchObj(yCurr)
            # println("Base Fun lineSearchObj -> $(size(lc)) -> $lc")
            lgc = lineSearchdfdx(yCurr)
            # println("Base Grad lineSearchdfdx -> $(size(lgc)) ->")
            # display(lgc)
            rlgc = dk'lgc
            # println("Results include $(size(rlgc)) ->")
            # display(rlgc)

            # Attempt a linesearch
            # Get the line search recommendation
            y0New, stepLS = backtrackLineSearch(yCurr, dk,
                            lineSearchObj, lineSearchdfdx, sp.paramA, sp.paramB)

            if true
                println("    Recommended Line Search Step: $stepLS")
                if stepLS < (sp.paramB ^ 3)
                    println("    Very low step.")
                end
            end

            if verbose
                print("Expected x = ")
                println("$y0New ?= $(primalVec(yCurr) + stepLS * dk)\n")
            end

            # Now we update the trust region
            yChange = norm(stepLS * dk)
            trustDelta = (yChange + sp.trc4 * trustDelta) / 2

        end

        # Save the new state to the output list
        push!(yNewtStates, y0New)

        if true
            println("Added State")#"$y0New\n")
        end

        # Break by tolerance and trust region size
        if (norm(residual, 2) < sp.rTol)
            println("Ended from tolerance at $i Newton steps\n")
            early = true
            break
        end

        # Update the current state with the new state
        yCurr = y0New

        if typeof(al.cM) == constraintManager_Dynamics
            aSize = size(al.cM.affineLambdaList, 1)
            ySize = size(y0New, 1)
            al.cM.affineLambdaList = y0New[ySize - aSize + 1:end]
        end

    end

    if !early
        println("Ended from max steps in $(sp.maxNewtonSteps) Newton Steps\n")
    end

    if verbose
        println("Ended Newton Method at $yCurr")
        println("ϕ(y) = $(lineSearchObj(primalVec(yCurr)))")
        println("∇ϕ(y) = $(lineSearchdfdx(primalVec(yCurr)))")
        cCurr = evalConstraints(al.cM, yCurr, al.rho)
        println("Constraints: $cCurr\n")
    end

    return yNewtStates, residNewt

end

"""
    ALPrimalNewtonMain(y0::SOCP_primals, al::augLagQP_2Cone,
                           sp::solverParams, verbose = false)

Performs a maximum on `N` outer steps as specified by [`solverParams`](@ref).
May return early if the residual is lower than the tolerance specified by
[`solverParams`](@ref) and the constraint penalty is high enough.

The function operates on an augmented lagrangian (AL) with primals (P) only
(as opposed to a Primal-Dual set-up) for second-order cone programs (SOCP).
See [`newtonTRLS_ALP`](@ref) for more information on the Newton Steps

returns (yStates, residuals) which are:
- `yStates`: An array with one entry of the primal vector per Newton step
- `residuals`: An array with one entry of the residual vector per Newton step

See also: [`augLagQP_2Cone`](@ref), [`SOCP_primals`](@ref)
"""
function ALPrimalNewtonMain(y0, al::augLag, sp::solverParams, verbose = false)

    yStates = []
    residuals = []
    push!(yStates, y0)

    for i in 1:(sp.maxOuterIters)

        # Update x at each iteration
        if verbose
            println("\n--------------------------------------")
            println("Next Full Update starting at $y0")
        end

        (yNewStates, resAtStates) = newtonTRLS_ALP(y0, al, sp, verbose)

        # Take each step in the arrays above and save it to the respective
        # overall arrays
        yStates = [yStates; yNewStates]
        residuals = [residuals; resAtStates]

        # Determine the new lambda and rho
        # λ ← λ + ρ c(x_k*)       - Which is to say update with prior x*
        # ρ ← min(ρ * 10, 10^6)   - Which is to say we bound ρ's growth by 10^6
        yNewest = yNewStates[end]

        yNewPrimals = getPrimals(al, yNewest)

        cCurr = evalConstraints(al.cM, yNewPrimals, al.rho)

        # Only turn on if the solver is not returning all of the points.
        # if false
        #     println()
        #     println("yStates = $yStates")
        #     println("yNewest = $yNewest")
        #     println("All residuals = $residuals")
        # end
        if verbose
            println("\n################")
            println("Former Lambda: $(al.lambda)")
        end

        updateDual!(al.cM, yNewPrimals, al.rho)
        al.rho = clamp(al.rho * sp.penaltyStep, 0, sp.penaltyMax)

        if verbose
            println("New state added: (Newest) $yNewest")
            println("cCurr: $(cCurr)")
            println("Lambda Updated: $(al.lambda) vs. $lambdaNew")
            println("rho updated: $(al.rho)")
            println("################\n")
        end

        if (norm(residuals[end], 2) < sp.rTol) && (al.rho == sp.penaltyMax)
            println("Ended early at $i outer steps")
            break
        else
            y0 = yNewest
        end
    end

    println("\nSOLVE COMPLETE\n")

    return yStates, residuals
end
