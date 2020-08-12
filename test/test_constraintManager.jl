# Now test the constraint on the constraint manager

using TrajOptSOCPs
using Test

@testset "Simple Constraint Manager Base" begin

    c1 = AL_AffineEquality([5], [10])
    c2 = AL_AffineEquality([5], [10])
    c3 = AL_AffineEquality([5], [10])
    lambdaVec = [[2], [3], [4]]

    cM1 = TrajOptSOCPs.constraintManager_Base([c1, c2, c3], lambdaVec)
    x0Test = [5]
    penaltyTest = 10.0

    # Evaluate the constraints
    @test evalConstraints(cM1, x0Test, penaltyTest) == 297.0

    # Evaluate Jacobian of the Constraint Terms
    gcm = evalGradConstraints(cM1, x0Test, penaltyTest)
    @test size(gcm) == ()

end
