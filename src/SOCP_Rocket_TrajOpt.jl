# This file will contain the SOCP TrajOpt Module
"""
    SOCP_Rocket_TrajOpt

A native Julia library to solve trajectory optimization problems that
contain second-order cone constraints.

Author: Daniel Neamati (Summer 2020)

Funding graciously provided by Caltech SURF program and the Homer J. Stewart
Fellowship.
"""
module SOCP_Rocket_TrajOpt

using LinearAlgebra, SparseArrays, Plots

# Declare Abstract Types
abstract type constraint end
abstract type objectiveFunc end

# Include all of the files
include("rocket\\rocket-setup.jl")
include("rocket\\ground.jl")
include("rocket\\maxThrust.jl")


include("dynamics\\trajectory-setup.jl")
include("objective\\LQR_objective.jl")
include("constraints\\constraintManager.jl")
include("constraints\\constraintManagerDynamics.jl")
include("auglag\\auglag-core.jl")

include("solver\\AL-Primal-Main-Solver.jl")
include("solver\\QP-AffineEquality-Solver.jl")
include("other_utils\\parsePrimalDual.jl")

include("results\\trajectoryParsing.jl")
include("results\\plotTrajectory.jl")
include("results\\plotConstraintViolation.jl")
include("results\\plotObjective.jl")
include("results\\batchPlots.jl")

# Export Structs
export augLag, solverParams, constraintManager_Dynamics, LQR_QP_Referenced,
       rocket_simple

# Export Functions
export ALPrimalNewtonMain, makeLQR_TrajReferenced

end  # module SOCP_TrajOpt
