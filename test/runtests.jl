using TrajOptSOCPs
using Test

@testset "Base Cases" begin
    # Write your tests here.
    @test fCheck2(4) == 16
    @test fCheck2(8) == 64
    @test fCheck2(10) == 100
end
