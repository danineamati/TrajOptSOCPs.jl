using SafeTestsets

# The absolute most basic tests
@safetestset "Most Basic Tests" begin
    include("test_mostBasic.jl")
end

# The affine constraint tests
@safetestset "Affine Constraints Tests" begin
    include("test_affineConstraints.jl")
end

# The base constraint manager tests
@safetestset "Base Constraint Manager Tests" begin
    include("test_constraintManager.jl")
end
  
 # Test the angle cone
@safetestset "Angle SOCP Constraint Test" begin
    include("test_angleSOCP.jl")
end
