include("../so3/rotation_matrix.jl")

using Random: rand, seed!
using Test: @testset, @test, @test_throws

seed!(0)

X90DEG_ARR = [1 0  0;
              0 0 -1;
              0 1  0]
Y90DEG_ARR = [0 0 1;
              0 1 0;
             -1 0 0]
Z90DEG_ARR = [0 -1 0;
              1  0 0;
              0  0 1]
X90DEG = RotationMatrix(X90DEG_ARR)
Y90DEG = RotationMatrix(Y90DEG_ARR)
Z90DEG = RotationMatrix(Z90DEG_ARR)

@testset "RotationMatrix Construction" begin
    # Test construction, with non-rotation matrix array but with checks turned off.
    @test isa(RotationMatrix(rand(2, 2), checks=false), RotationMatrix)
    @test isa(RotationMatrix(rand(3, 3), checks=false), RotationMatrix)
    @test isa(RotationMatrix(rand(4, 4), checks=false), RotationMatrix)

    # Test construction, with non-rotation matrix array with checks turned on.
    @test_throws ArgumentError RotationMatrix(rand(2, 2))
    @test_throws ArgumentError RotationMatrix(rand(3, 3))
    @test_throws ArgumentError RotationMatrix(rand(4, 4))

    # Test construction with proper rotation matrices
    @test isa(RotationMatrix([1 0 0; 0 1 0; 0 0 1]), RotationMatrix)
    @test isa(RotationMatrix(X90DEG_ARR), RotationMatrix)
    @test isa(RotationMatrix(Y90DEG_ARR), RotationMatrix)
    @test isa(RotationMatrix(Z90DEG_ARR), RotationMatrix)

end

@testset "RotationMatrix Identity" begin
    eye = one(RotationMatrix)
    @test eye * eye ≈ eye
    @test eye * X90DEG ≈ X90DEG
    @test eye * Y90DEG ≈ Y90DEG
    @test eye * Z90DEG ≈ Z90DEG
    @test X90DEG * eye ≈ X90DEG
    @test Y90DEG * eye ≈ Y90DEG
    @test Z90DEG * eye ≈ Z90DEG
end

@testset "RotationMatrix Multiplication" begin
    # RotationMatrix-RotationMatrix multiplication
    @test X90DEG * Y90DEG ≈ RotationMatrix([0 0 1; 1 0 0; 0 1 0])
    @test X90DEG * Z90DEG ≈ RotationMatrix([0 -1 0; 0 0 -1; 1 0 0])
    @test Y90DEG * X90DEG ≈ RotationMatrix([0 1 0; 0 0 -1; -1 0 0])
    @test Y90DEG * Z90DEG ≈ RotationMatrix([0 0 1; 1 0 0; 0 1 0])
    @test Z90DEG * X90DEG ≈ RotationMatrix([0 0 1; 1 0 0; 0 1 0])
    @test Z90DEG * Y90DEG ≈ RotationMatrix([0 -1 0; 0 0 1; -1 0 0])

    # RotationMatrix-vector multiplication
    v = [1; 2; 3]
    @test X90DEG * v ≈ [1; -3; 2]
    @test Y90DEG * v ≈ [3; 2; -1]
    @test Z90DEG * v ≈ [-2; 1; 3]

    # RotationMatrix-matrix multiplication
    mat = [1 2 3;
           2 3 1;
           3 1 2]
    @test X90DEG * mat ≈ [1  2  3;
                         -3 -1 -2;
                          2  3  1]
    @test Y90DEG * mat ≈ [3  1  2;
                          2  3  1;
                         -1 -2 -3]
    @test Z90DEG * mat ≈ [-2 -3 -1;
                           1  2  3;
                           3  1  2]
end

@testset "RotationMatrix Inverse" begin
    eye = one(RotationMatrix)
    @test inv(eye) ≈ eye
    @test eye * inv(eye) ≈ eye
    @test inv(eye) * eye ≈ eye

    @test inv(X90DEG) * X90DEG ≈ eye 
    @test inv(Y90DEG) * Y90DEG ≈ eye 
    @test inv(Z90DEG) * Z90DEG ≈ eye 
    @test X90DEG * inv(X90DEG) ≈ eye 
    @test Y90DEG * inv(Y90DEG) ≈ eye 
    @test Z90DEG * inv(Z90DEG) ≈ eye 

end

@testset "RotationMatrix Logarithmic Map" begin
    @test log(one(RotationMatrix)) ≈ zero(SO3Algebra)
    @test log(X90DEG) ≈ SO3Algebra([0.5 * π; 0; 0])
    @test log(Y90DEG) ≈ SO3Algebra([0; 0.5 * π; 0])
    @test log(Z90DEG) ≈ SO3Algebra([0; 0; 0.5 * π])
end

@testset "RotationMatrix Conversion" begin
    eye = one(RotationMatrix)

    @test eye ≈ convert(RotationMatrix, eye)
    @test X90DEG ≈ convert(RotationMatrix, X90DEG)
    @test Y90DEG ≈ convert(RotationMatrix, Y90DEG)
    @test Z90DEG ≈ convert(RotationMatrix, Z90DEG)

    test_mat = RotationMatrix(rand(2, 2), checks=false)
    @test_throws ArgumentError convert(RotationMatrix, test_mat, checks=true)
    @test test_mat ≈ convert(RotationMatrix, test_mat, checks=false)

    @test convert(Array, eye) ≈ eye.value
end

@testset "RotationMatrix Correctness Checks" begin
    eye = one(RotationMatrix)

    @test check(eye)[1]
    @test check(convert(Array, eye))[1]
    @test !check(RotationMatrix(rand(3, 3), checks=false))[1]
    @test !check(rand(3, 3))[1]
end

@testset "RotationMatrix Approximation Check" begin
    eye = one(RotationMatrix)

    @test eye ≈ eye
    @test eye ≉ RotationMatrix(rand(3, 3), checks=false)
end
