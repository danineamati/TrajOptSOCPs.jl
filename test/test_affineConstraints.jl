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
end

@testset "Affine Equality Constraint" begin
    # Create a simple AL_AffineEquality object
    dT2 = AL_AffineEquality([5 0; 0 6], [10; 12])

    # Check that it can determine when it is on the line
    @test TrajOptSOCPs.satisfied(dT2, [2; 2])

    # And likewise when it is not on the line
    ws = TrajOptSOCPs.whichSatisfied(dT2, [2; 4]) #(should return (True, False))

    @test ws[1]
    @test !ws[2]

    # Test Projections and Violations
    testpt = [10; 12]
    # This point is not on the intersetion:
    @test !TrajOptSOCPs.satisfied(dT2, testpt)

    # Get the projections and norm
    # pv = getProjVecs(dT2, testpt)
    gN = TrajOptSOCPs.getNormToProjVals(dT2, testpt)

    @test (gN'gN) != 0.0

    # Test the Calculus
    gradDT2 = TrajOptSOCPs.getGradC(dT2, testpt)
    hessDT2 = TrajOptSOCPs.getHessC(dT2)

    # Here we simply check the size.
    @test size(gradDT2, 1) == 2
    @test size(hessDT2) == (2, 2)

end

    # println("\nAffine Inequality Constraints")
    # diT = AL_AffineInequality([5], [10])
    # println("** Simplest")
    # println("On line    (True) : $(TrajOptSOCPs.satisfied(diT, 2))")
    # println("Above line (False): $(satisfied(diT, 4))")
    # println("Under line (True) : $(satisfied(diT, 0))")
    # diT2 = AL_AffineInequality([5 0; 0 6], [10; 12])
    # println("** 2D")
    # println("On Intesection (True)        : $(satisfied(diT2, [2; 2]))")
    # println("Check which not (True, False): $(whichSatisfied(diT2, [2; 4]))")
    # println("Check which not (True, True) : $(whichSatisfied(diT2, [2; 0]))")
    # println("** Testing the projections and violations")
    # testpt = [10; 12]
    # println("Point $testpt not at in region : $(!satisfied(diT2, testpt))")
    # println("Projection onto feasible region : $(getProjVecs(diT2, testpt))")
    # gN = getNormToProjVals(diT2, testpt)
    # println("Violation is : $gN")
    # println("Total Violation is : $(gN'gN)")
    # testpt2 = [-10; -12]
    # println("Point $testpt2 in feasible region : $(satisfied(diT2, testpt2))")
    # println("Projection onto feasible region : $(getProjVecs(diT2, testpt2))")
    # gN2 = getNormToProjVals(diT2, testpt2)
    # println("Violation is : $gN2")
    # println("Total Violation is : $(gN2'gN2)")
    # println("** Testing Calculus")
    # println("Hessian = $(getHessC(diT2))")
    # println("Grad at $testpt = $(getGradC(diT2, testpt))")
    # println("Grad at $testpt2 = $(getGradC(diT2, testpt2))")
