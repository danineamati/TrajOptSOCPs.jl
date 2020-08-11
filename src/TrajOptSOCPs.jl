# This file will contain the SOCP TrajOpt Module
"""
    TrajOptSOCPs

A native Julia library to solve trajectory optimization problems that
contain second-order cone constraints.

Author: Daniel Neamati (Summer 2020)

Funding graciously provided by Caltech SURF program and the Homer J. Stewart
Fellowship.
"""
module TrajOptSOCPs

using LinearAlgebra, SparseArrays

# Declare Abstract Types
abstract type constraint end
abstract type objectiveFunc end
abstract type constraintManager end

# Basic Test
include("extra.jl")
export fCheck2

# Include all of the files
include("rocket/rocket-setup.jl")
include("rocket/ground.jl")
include("rocket/maxThrust.jl")


include("dynamics/trajectory-setup.jl")
include("objective/LQR_objective.jl")
include("constraints/constraintManager.jl")
include("constraints/constraintManagerDynamics.jl")
include("auglag/auglag-core.jl")

include("solver/AL-Primal-Main-Solver.jl")
include("solver/QP-AffineEquality-Solver.jl")
include("other_utils/parsePrimalDual.jl")

include("results/trajectoryParsing.jl")
# include("results/plotTrajectory.jl")
# include("results/plotConstraintViolation.jl")
# include("results/plotObjective.jl")
# include("results/batchPlots.jl")

# # Export Constraint Structs
export AL_AffineEquality, AL_AffineInequality,
        AL_simpleCone, AL_Multiple_simpleCone

# Export Top Level Structs
export augLag, solverParams, constraintManager_Dynamics,
        LQR_QP_Referenced, rocket_simple

# Export Top Level Functions
export ALPrimalNewtonMain, solParamPrint,
        initializeTraj, rocketDynamicsFull

# Export Helper Set-Up Functions
export makeLQR_TrajReferenced, makeMaxThrustConstraint, makeGroundConstraint

# Export Evaluation functions
export evalConstraints, evalGradConstraints, evalHessConstraints,
       evalAL, evalGradAL, evalHessAl

# Export Helper Utils
export parsePrimalDualVec

end  # module SOCP_TrajOpt
