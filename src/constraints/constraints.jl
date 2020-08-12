
using SparseArrays


abstract type constraint end

# Constraint Types and settings

include("affineConstraints.jl")
include("socpConstraints.jl")
include("socpNewConstraints.jl")
include("socpAngleConstraints.jl")


#= Augmented Lagrangian Constraints
These are presented in the order:
- AffineEquality
- AffineInequality
- (Old) SOCP Constraints
- (Old) Functional SOCP Constraints with Slack Variables
- Simple SOCP Constraints
- Simple SOCP Constraints for N-constraints
- SOCP Constraints for Angles
- SOCP Constraints for Angles for N-constraints
=#
