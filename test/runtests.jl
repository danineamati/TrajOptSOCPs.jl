using SafeTestsets

# The absolute most basic tests
@safetestset "Most Basic Tests" begin
    include("test_mostBasic.jl")
end

# The affine constraint tests
@safetestset "Affine Constraints Tests" begin
    include("test_affineConstraints.jl")
end
