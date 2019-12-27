include("../so3/so3.jl")

using Test: @testset

@testset "RotationMatrix" begin
    include("test_rotation_matrix.jl")
end

@testset "SO3Algebra" begin
    include("test_so3_algebra.jl")
end
