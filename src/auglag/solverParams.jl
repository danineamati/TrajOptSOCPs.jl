# Solver Parameters

# -------------------------
# Solver Parameters
# -------------------------
"""
Line Search Parameters (Notation from Toussaint Notes (2020))

`paramA::Float16`         # Used in Line Search, should be [0.01, 0.3]

`paramB::Float16`         # Used in Line Search, should be [0.1, 0.8]

----
Iteration Counter Parameters

`maxOuterIters::Int32`    # Number of Outer Loop iterations

`maxNewtonSteps::Int32`   # Number of Newton Steps per Outer Loop iterations

----
Exit Condition Parameter

`rTol::Float64`           # When steps are within rTol, loop will stop.

----
Penalty Update Parameters

`penaltyStep::Float16`    # Multiplies the penalty parameter per outer loop

`penaltyMax::Float64`     # Maximum value of the penalty parameter

----
Trust Region Parameters (Notation from Nocedal et Yuan (1998))

`trSizeStart::Float32`    # Starting size of the trust region

`trc1::Float16`           # Success Size Increase Parameter (1 < c1)

`trc2::Float16`           # Taylor Series Error Parameter (0 < c2 < 1)

`trc3::Float16`           # Failed Size Reduction Parameter (0 < c3 < c4)

`trc4::Float16`           # Failed Size Reduction Parameter (c3 < c4 < 1)
"""
struct solverParams
    # Line Search Parameters (Notation from Toussaint Notes (2020))
    paramA::Float16         # Error allowed should be [0.01, 0.3]
    paramB::Float16         # Reduction Factor should be [0.1, 0.8]
    # Iteration Counter Parameters
    maxOuterIters::Int32    # Number of Outer Loop iterations
    maxNewtonSteps::Int32   # Number of Newton Steps per Outer Loop iterations
    # Exit Condition Parameter
    rTol::Float64           # When steps are within rTol, loop will stop.
    # Penalty Update Parameters
    penaltyStep::Float16    # Multiplies the penalty parameter per outer loop
    penaltyMax::Float64     # Maximum value of the penalty parameter
    # Trust Region Parameters (Notation from Nocedal et Yuan (1998))
    trSizeStart::Float32    # Starting size of the trust region
    trc1::Float16           # Success Size Increase Parameter (1 < c1)
    trc2::Float16           # Taylor Series Error Parameter (0 < c2 < 1)
    trc3::Float16           # Failed Size Reduction Parameter (0 < c3 < c4)
    trc4::Float16           # Failed Size Reduction Parameter (c3 < c4 < 1)
end

"""
Prints the struct variables in the `solverParams` struct.

See also: [`solverParams`](@ref) for more information.
"""
function solParamPrint(sp::solverParams)
    println()
    println("Beginning solver with parameters: ")
    println("(Line Search) : a = $(sp.paramA), b = $(sp.paramB)")
    print("(Loop #)      : Outer = $(sp.maxOuterIters), ")
    println("Newton = $(sp.maxNewtonSteps)")
    println("(or End at)   : Δ(∇L) = $(sp.rTol)")
    println("(Penalty)     : Δρ = $(sp.penaltyStep), ρMax = $(sp.penaltyMax)")
    println("(Trust Region): Δ = $(sp.trSizeStart)")
    println("                c1 = $(sp.trc1) (with 1 < c1 -> $(1 < sp.trc1))")
    print("                c2 = $(sp.trc2) (with 0 < c2 < 1 -> ")
    println("$(0 < sp.trc2) & $(sp.trc2 < 1))")
    println("                c3 = $(sp.trc3), c4 = $(sp.trc4) ")
    print("                (with 0 < c3 < c4 < 1) -> $(0 < sp.trc3) & ")
    println("$(sp.trc3 < sp.trc4) & $(sp.trc4 < 1)")
end
