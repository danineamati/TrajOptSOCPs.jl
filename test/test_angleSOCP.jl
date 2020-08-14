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

@testset "Test Values of Angle SOCP" begin
    sI = Diagonal([1; 0])
    tI = [0; 1]

    alpha = TrajOptSOCPs.getAlpha(20, true)
    @test alpha == -tand(20)

    asocp = AL_simpleAngleCone(alpha, sI, tI)

    xTest1 = [1; -1] # Angle of 45 deg > 20 deg
    xTest2 = [0.18; -1] # Angle of ~10 deg < 20 deg
    xTest3 = [alpha; -1] # Exactly 20 deg

    @test !TrajOptSOCPs.satisfied(asocp, xTest1)
    @test TrajOptSOCPs.satisfied(asocp, xTest2)
    @test TrajOptSOCPs.satisfied(asocp, xTest3)

    # Test the constraint violation
    gN1 = TrajOptSOCPs.getNormToProjVals(asocp, xTest1)
    gN2 = TrajOptSOCPs.getNormToProjVals(asocp, xTest2)
    gN3 = TrajOptSOCPs.getNormToProjVals(asocp, xTest3)

    @test gN1'gN1 != 0.0
    @test gN2'gN2 == 0.0
    @test gN3'gN3 == 0.0

end

@testset "Test Calculus of Angle SOCP" begin
    sI = Diagonal([1; 0])
    tI = [0; 1]

    alpha = TrajOptSOCPs.getAlpha(20, true)

    asocp = AL_simpleAngleCone(alpha, sI, tI)

    xTest1 = [1; -1] # Angle of 45 deg > 20 deg
    xTest2 = [0.18; -1] # Angle of ~10 deg < 20 deg
    xTest3 = [alpha; -1] # Exactly 20 deg

    # Test the constraint violation
    g1 = TrajOptSOCPs.getGradC(asocp, xTest1)
    g2 = TrajOptSOCPs.getGradC(asocp, xTest2)
    g3 = TrajOptSOCPs.getGradC(asocp, xTest3)

    @test size(g1, 1) == 2
    @test size(g2, 1) == 2
    @test size(g3, 1) == 2

    @test g1 == [1; -alpha]
    @test g2 == [0.0; 0.0] #[1; -alpha]
    @test g3 == [0.0; 0.0] #[sign(alpha); -alpha]

    h1 = TrajOptSOCPs.getHessC(asocp, xTest1)
    h2 = TrajOptSOCPs.getHessC(asocp, xTest2)
    h3 = TrajOptSOCPs.getHessC(asocp, xTest3)

    @test h1 == h2
    @test h2 == h3
    @test h3 == spzeros(2, 2)

    ah1 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest1)
    ah2 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest2)
    ah3 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest3)

    @test ah1 == [1 -alpha; -alpha alpha^2]
    @test ah2 == ah3
    @test h3 == ah3

end

@testset "Test Calculus of Angle SOCP 3D" begin
    sI = Diagonal([1; 0; 1])
    tI = [0; 1; 0]

    alpha = TrajOptSOCPs.getAlpha(20, true)

    asocp = AL_simpleAngleCone(alpha, sI, tI)

    xTest1 = [1; -1; 1] # Angle of 45 deg > 20 deg
    xTest2 = [0.125; -1; 0.125] # Angle of ~10 deg < 20 deg
    xTest3 = [0; -1; alpha] # Exactly 20 deg
    xTest4 = [0.1; -0.2; 0.3] # Angle of ~58 deg > 20 deg

    @test !TrajOptSOCPs.satisfied(asocp, xTest4)

    # Test the constraint violation
    g1 = TrajOptSOCPs.getGradC(asocp, xTest1)
    g2 = TrajOptSOCPs.getGradC(asocp, xTest2)
    g3 = TrajOptSOCPs.getGradC(asocp, xTest3)
    g4 = TrajOptSOCPs.getGradC(asocp, xTest4)

    @test size(g1, 1) == 3
    @test size(g2, 1) == 3
    @test size(g3, 1) == 3
    @test size(g4, 1) == 3

    @test g1 == [1 / sqrt(2); -alpha; 1 / sqrt(2)]
    @test g2 == zeros(3) #[1 / sqrt(2); -alpha; 1 / sqrt(2)]
    @test g3 == zeros(3) #[0; -alpha; sign(alpha)]
    @test g4 == [0.1 / norm([0.1; 0.3]); -alpha; 0.3 / norm([0.1; 0.3])]

    h1 = TrajOptSOCPs.getHessC(asocp, xTest1)
    h2 = TrajOptSOCPs.getHessC(asocp, xTest2)
    h3 = TrajOptSOCPs.getHessC(asocp, xTest3)
    h4 = TrajOptSOCPs.getHessC(asocp, xTest4)

    @test h1 == h2
    @test h2 == h3
    @test h3 == h4
    @test h4 == spzeros(3, 3)

    ah1 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest1)
    ah2 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest2)
    ah3 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest3)
    ah4 = TrajOptSOCPs.getHessC_ALTerm(asocp, xTest4)

    sq2 = - alpha / sqrt(2)
    ah1_check = [1/2 sq2 1/2; sq2 alpha^2 sq2; 1/2 sq2 1/2]

    ns4 = norm([0.1; 0.3])
    sq4 = - alpha / norm([0.1; 0.3])

    ah4_check = [((0.1^2) / ns4^2) (0.1 * sq4) ((0.1 * 0.3) / ns4^2);
                (0.1 * sq4)         alpha^2         (0.3 * sq4);
                ((0.1 * 0.3) / ns4^2) (0.3 * sq4) ((0.3^2) / ns4^2)]

    @test round.(ah1, digits = 8) == round.(ah1_check, digits = 8)
    @test ah2 == ah3
    @test h3 == ah3
    @test round.(ah4, digits = 8) == round.(ah4_check, digits = 8)

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

    @test size(gN, 1) == 3
    @test gN'gN != 0.0

    jacobMany = TrajOptSOCPs.getGradC(aMany, xTest)
    hessMany = TrajOptSOCPs.getHessC(aMany, xTest)
    ahessMany = TrajOptSOCPs.getHessC_ALTerm(aMany, xTest)

    xSize = size(xTest, 1)

    @test size(jacobMany) == (3, xSize)
    @test size(hessMany) == (xSize, xSize)
    @test size(ahessMany) == (xSize, xSize)

end
