using TrajectoryOptimizationWithSOCPs
using Test

@testset "TrajectoryOptimizationWithSOCPs.jl" begin
    # Write your tests here.
    @test fCheck2(4) == 16
    @test fCheck2(8) == 64
end
