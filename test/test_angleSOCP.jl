# Tests for angle SOCP constriants

using TrajOptSOCPs
using Test
using SparseArrays, LinearAlgebra

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


@testset "Multiple Angle SOCP" begin
    sI = Diagonal([1; 0])
    tI = [0; 1]
    alpha = 0.2

    nDim = 2
    indicator = findnz(sparse(
                    [0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 1; 0; 0; 0; 0; 0; 1; 0]))[1]

    aMany = TrajOptSOCPs.makeAL_Multiple_AngleCone(alpha, sI, tI,
                                                    nDim, indicator)

    # Check that it worked
    @test typeof(aMany) == AL_Multiple_AngleCone

    xTest = [1; 2; 3; 4; 5; 6; 6; 5; 4; 3; 2; 1; 1; 2; 3; 4; 5; 6]

    # The code needs to notice that there are three vectors to extract and then
    # get those and evaluate it.
    rawVec = TrajOptSOCPs.getRaw(aMany, xTest)

    @test size(rawVec, 1) == 3

    gN = TrajOptSOCPs.getNormToProjVals(aMany, xTest)

    @test size(gN, 1) ==3
    @test gN'gN != 0.0

end
