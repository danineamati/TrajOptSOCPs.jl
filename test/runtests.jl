using SafeTestsets

# The absolute most basic tests
@safetestset "Most Basic Tests" begin
    include("test_mostBasic.jl")
end

# Test the angle cone
@safetestset "Angle SOCP Constraint Test" begin
    include("test_angleSOCP.jl")
end
