# Tests for angle SOCP constriants

using TrajOptSOCPs
using Test, LinearAlgebra

@testset "Simple Angle SOCP" begin
    sI = Diagonal([1; 0])
    tI = [0; 1]
    asocp = AL_simpleAngleCone(0.2, sI, tI)

    # The first and third should be satisified. The second should not be
    # satisfied.
    xTestArr = [[0.1; 2], [2; 2], [-0.4; 3]]
    ws = TrajOptSOCPs.whichSatisfied(asocp, xTestArr)

    @test ws[1]
    @test !ws[2]
    @test ws[3]

    gN1 = TrajOptSOCPs.getNormToProjVals(asocp, xTestArr[1])
    gN2 = TrajOptSOCPs.getNormToProjVals(asocp, xTestArr[2])

    @test gN1'gN1 == 0
    @test gN2'gN2 != 0


    # Now the calculus
    grad1 = TrajOptSOCPs.getGradC(asocp, xTestArr[1])
    hess1 = TrajOptSOCPs.getHessC(asocp, xTestArr[1])
    alhess1 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTestArr[1])

    @test size(grad1, 1) == 2
    @test size(hess1) == (2, 2)
    @test size(alhess1) == (2, 2)

    grad2 = TrajOptSOCPs.getGradC(asocp, xTestArr[2])
    hess2 = TrajOptSOCPs.getHessC(asocp, xTestArr[2])
    alhess2 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTestArr[2])

    @test size(grad2, 1) == 2
    @test size(hess2) == (2, 2)
    @test size(alhess2) == (2, 2)

end
