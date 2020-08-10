# Dynamically aware constraints manager

include("constraintManager.jl")
include("..\\objective\\LQR_objective.jl")


"""
    constraintManager_Dynamics

Contains a list of constraints and the dual variables (lambda) that correspond
with each constraint.

Unlike [`constraintManager_Base`](@ref), this constraint manager also has the
dynamics constraints embedded such that the nearest dynamically feasible
trajectory can be determined at any given time.
"""
mutable struct constraintManager_Dynamics <: constraintManager
    cList::Array{constraint, 1}
    lambdaList
    dynamicsConstraint::AL_AffineEquality
    affineLambdaList
end


#############################################################
###            Evaluate the Affine Constraints            ###
#############################################################

function evalAffineEq(cM::constraintManager_Dynamics, yCurr)
    return getRaw(cM.dynamicsConstraint, yCurr)
end

function affineBlock(cM::constraintManager_Dynamics)
    return cM.dynamicsConstraint.A
end
