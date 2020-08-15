# Test MonteCarlo functions

using TrajOptSOCPs
using Test

@testset "Feasible" begin

    # Define the polyhedron
    A = [4 5; -4 5; 0 -1]
    b = [20; 30; 1]

    # Test a point in the feasible region (true)
    @test TrajOptSOCPs.isFeasiblePolyHedron(A, b, [-1; 0.5])

    # Test three points outside the feasible region (false)
    @test !TrajOptSOCPs.isFeasiblePolyHedron(A, b, [-10; 0.5])
    @test !TrajOptSOCPs.isFeasiblePolyHedron(A, b, [-1; 10])
    @test !TrajOptSOCPs.isFeasiblePolyHedron(A, b, [-1; -10])
end

@testset "MonteCarlo" begin

    # Generating a 3D point
    pt1 = TrajOptSOCPs.getPoint([3; 2; -5], [10; 100; 0.5])

    # Checking that the point is in the bounds
    @test (10 > pt1[1] > 3)
    @test (100 > pt1[2] > 2)
    @test (0.5 > pt1[3] > -5)

    # Generating 3 Random 4D Points
    numpts = 3
    pts = TrajOptSOCPs.montecarlo([3; 2; -5; -0.01],
                                  [10; 100; 0.5; 0.01],
                                  numpts)

    @test size(pts, 1) == numpts

    # Check that the points are indeed 4D
    @test size(pts[1], 1) == 4
    @test size(pts[numpts], 1) == 4

    # Generating 2 Random Points in a smaller [0, 1] box:
    pts2 = TrajOptSOCPs.montecarlo([-20; -20], [20; 20], 2, TrajOptSOCPs.inBox)

    @test size(pts2, 1) == 2

    # Check that the points are indeed 2D
    @test size(pts2[1], 1) == 2
    @test size(pts2[2], 1) == 2

end
