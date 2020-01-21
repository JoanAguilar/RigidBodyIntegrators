include("../so3/so3.jl")
include("../se3/se3.jl")

using Test: @testset

@testset "RotationMatrix" begin
    include("test_rotation_matrix.jl")
end

@testset "SO3Algebra" begin
    include("test_so3_algebra.jl")
end

@testset "CompositeSE3Group" begin
    include("test_composite_se3_group.jl")
end
