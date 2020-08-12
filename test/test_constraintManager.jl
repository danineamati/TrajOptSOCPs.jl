# Now test the constraint on the constraint manager

using TrajOptSOCPs
using Test

@testset "Simple Constraint Manager Base" begin

    # Make three constraints...
    c1 = AL_AffineEquality([5], [10])
    c2 = AL_AffineEquality([5], [10])
    c3 = AL_AffineEquality([5], [10])
    # ... and the associated dual
    lambdaVec = [[2], [3], [4]]

    # Make the constraint manager
    cM1 = TrajOptSOCPs.constraintManager_Base([c1, c2, c3], lambdaVec)
    x0Test = [5]
    penaltyTest = 10.0

    # Evaluate the constraints
    @test evalConstraints(cM1, x0Test, penaltyTest) == 297.0

    # Evaluate Jacobian of the Constraint Terms
    gcm = evalGradConstraints(cM1, x0Test, penaltyTest)
    @test size(gcm) == ()

end


@testset "2D Equalities and Inequalities Mixed Case" begin

    # Make three constraints of mixed types...
    c1 = AL_AffineInequality([5 0; 0 6], [10; 12])
    c2 = AL_AffineEquality([5 0; 0 6], [10; 12])
    c3 = AL_AffineInequality([5 0; 0 6], [10; 12])
    # ... and the associated dual
    lambdaVec = [[2; 2], [3; 4], [4; 2]]

    # Make the constraint manager
    cM1 = TrajOptSOCPs.constraintManager_Base([c1, c2, c3], lambdaVec)
    x0Test = [5; -3]
    penaltyTest = 10.0

    # Evaluate the constraints
    @test evalConstraints(cM1, x0Test, penaltyTest) == 567.0

    # Evaluate Gradient of the Constraint Terms
    gcm = evalGradConstraints(cM1, x0Test, penaltyTest)
    @test size(gcm, 1) == 2


    # Evaluate Hessian of the Constraint Terms
    hcm = evalHessConstraints(cM1, x0Test)
    @test size(hcm) == (2, 2)

end
