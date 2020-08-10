using TrajectoryOptimizationWithSOCPs
using Test

@testset "TrajectoryOptimizationWithSOCPs.jl" begin
    # Write your tests here.
    @test fCheck(4) == 16
    @test fCheck(8) == 64
end
