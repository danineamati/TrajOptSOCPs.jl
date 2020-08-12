# Test Affine Constraints

using TrajOptSOCPs
using Test

@testset "Affine Equality Constraint Simple" begin
    # Create a simple AL_AffineEquality object
    dT = AL_AffineEquality([5], [10])

    # Check that it can determine when it is on the line
    @test TrajOptSOCPs.satisfied(dT, 2)

    # And likewise when it is not on the line
    @test !TrajOptSOCPs.satisfied(dT, 4)
    @test !TrajOptSOCPs.satisfied(dT, 0)

    # Run a simple evaluation
    @test TrajOptSOCPs.getRaw(dT, 2) == [0]
end

@testset "Affine Equality Constraint" begin
    # Create a higher dimensional AL_AffineEquality object
    dT2 = AL_AffineEquality([5 0; 0 6], [10; 12])

    # Check that it can determine when it is on the intersection
    @test TrajOptSOCPs.satisfied(dT2, [2; 2])

    # And likewise when it is not on the intersection
    ws = TrajOptSOCPs.whichSatisfied(dT2, [2; 4])
    #(should return (True, False))

    @test ws[1]
    @test !ws[2]

    # Test Projections and Violations
    testpt = [10; 12]
    # This point is not on the intersection:
    @test !TrajOptSOCPs.satisfied(dT2, testpt)

    # Get the projections and norm
    pv = TrajOptSOCPs.getProjVecs(dT2, testpt)
    gN = TrajOptSOCPs.getNormToProjVals(dT2, testpt)

    @test (gN'gN) != 0.0

    # Test the Calculus
    gradDT2 = TrajOptSOCPs.getGradC(dT2, testpt)
    hessDT2 = TrajOptSOCPs.getHessC(dT2)

    # Here we simply check the size.
    @test size(gradDT2, 1) == 2
    @test size(hessDT2) == (2, 2)

    # Test the the augmented lagrangian term also works
    alHess = TrajOptSOCPs.getHessC_ALTerm(dT2, testpt)

    @test size(alHess) == (2, 2)

end

@testset "Affine Inequality Constraint Simple" begin
    # Create a simple AL_AffineEquality object
    diT = AL_AffineInequality([5], [10])

    # Check that it can determine when it is on the line
    @test TrajOptSOCPs.satisfied(diT, 2)

    # And likewise when it is not on the line.
    # Above the line (false)
    @test !TrajOptSOCPs.satisfied(diT, 4)
    # Under the line (true)
    @test TrajOptSOCPs.satisfied(diT, 0)

    # Run a simple evaluation
    @test TrajOptSOCPs.getRaw(diT, 2) == [0]
end

@testset "Affine Inequality Constraint" begin
    # Create a higher dimensional AL_AffineInequality object
    diT2 = AL_AffineInequality([5 0; 0 6], [10; 12])

    # Check that it can determine when it is on the intersection
    @test TrajOptSOCPs.satisfied(diT2, [2; 2])

    # And likewise when it is not on the intersection
    ws = TrajOptSOCPs.whichSatisfied(diT2, [2; 4])
    #(should return (True, False))

    @test ws[1]
    @test !ws[2]

    # try a different set
    ws = TrajOptSOCPs.whichSatisfied(diT2, [2; 0])

    @test ws[1]
    @test ws[2]

    # Test Projections and Violations
    testpt = [10; 12]
    # This point is not on the intersection:
    @test !TrajOptSOCPs.satisfied(diT2, testpt)

    # Get the projections and violation
    pv = TrajOptSOCPs.getProjVecs(diT2, testpt)
    gN = TrajOptSOCPs.getNormToProjVals(diT2, testpt)

    @test (gN'gN) != 0.0


    # Try another test point
    testpt2 = [-10; -12]

    # This one should satisfy the constraint
    @test TrajOptSOCPs.satisfied(diT2, testpt2)

    # Get the projections and violation
    pv2 = TrajOptSOCPs.getProjVecs(diT2, testpt2)
    gN2 = TrajOptSOCPs.getNormToProjVals(diT2, testpt2)

    @test (gN2'gN2) == 0.0

    # Test the Calculus
    gradDT2 = TrajOptSOCPs.getGradC(diT2, testpt)
    gradDT22 = TrajOptSOCPs.getGradC(diT2, testpt2)
    hessDT2 = TrajOptSOCPs.getHessC(diT2)

    # Here we simply check the size.
    @test size(gradDT2, 1) == 2
    @test size(gradDT22, 1) == 2
    @test size(hessDT2) == (2, 2)

    # Test the the augmented lagrangian term also works
    alHess = TrajOptSOCPs.getHessC_ALTerm(diT2, testpt)

    @test size(alHess) == (2, 2)

end
