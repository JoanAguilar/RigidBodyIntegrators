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
EYE = one(RotationMatrix)
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
    @test EYE * EYE ≈ EYE
    @test EYE * X90DEG ≈ X90DEG
    @test EYE * Y90DEG ≈ Y90DEG
    @test EYE * Z90DEG ≈ Z90DEG
    @test X90DEG * EYE ≈ X90DEG
    @test Y90DEG * EYE ≈ Y90DEG
    @test Z90DEG * EYE ≈ Z90DEG
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
    @test inv(EYE) ≈ EYE
    @test EYE * inv(EYE) ≈ EYE
    @test inv(EYE) * EYE ≈ EYE

    @test inv(X90DEG) * X90DEG ≈ EYE
    @test inv(Y90DEG) * Y90DEG ≈ EYE
    @test inv(Z90DEG) * Z90DEG ≈ EYE
    @test X90DEG * inv(X90DEG) ≈ EYE
    @test Y90DEG * inv(Y90DEG) ≈ EYE
    @test Z90DEG * inv(Z90DEG) ≈ EYE
end

@testset "RotationMatrix Logarithmic Map" begin
    @test log(one(RotationMatrix)) ≈ zero(SO3Algebra)
    @test log(X90DEG) ≈ SO3Algebra([0.5 * π; 0; 0])
    @test log(Y90DEG) ≈ SO3Algebra([0; 0.5 * π; 0])
    @test log(Z90DEG) ≈ SO3Algebra([0; 0; 0.5 * π])
end

@testset "RotationMatrix Conversion" begin
    @test EYE ≈ convert(RotationMatrix, EYE)
    @test X90DEG ≈ convert(RotationMatrix, X90DEG)
    @test Y90DEG ≈ convert(RotationMatrix, Y90DEG)
    @test Z90DEG ≈ convert(RotationMatrix, Z90DEG)

    test_mat = RotationMatrix(rand(2, 2), checks=false)
    @test_throws ArgumentError convert(RotationMatrix, test_mat, checks=true)
    @test test_mat ≈ convert(RotationMatrix, test_mat, checks=false)

    @test convert(Array, EYE) ≈ EYE.value
end

@testset "RotationMatrix Correctness Checks" begin
    @test check(EYE)[1]
    @test check(convert(Array, EYE))[1]
    @test !check(RotationMatrix(rand(3, 3), checks=false))[1]
    @test !check(rand(3, 3))[1]
end

@testset "RotationMatrix Approximation Check" begin
    @test EYE ≈ EYE
    @test EYE ≉ RotationMatrix(rand(3, 3), checks=false)
end
