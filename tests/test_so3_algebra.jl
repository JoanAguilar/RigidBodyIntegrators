using Test: @testset, @test, @test_throws

# Instantiation of the zero element.
OH = zero(SO3Algebra)

# Unit elements along each coordinate, used in some tests.
ALG_X = SO3Algebra([1; 0; 0])
ALG_Y = SO3Algebra([0; 1; 0])
ALG_Z = SO3Algebra([0; 0; 1])

@testset "SO3Algebra Construction" begin
    # Test construction, with vectors that are not 3-dimensional but with checks turned off.
    @test isa(SO3Algebra(rand(2), checks=false), SO3Algebra)
    @test isa(SO3Algebra(rand(4), checks=false), SO3Algebra)
    @test isa(SO3Algebra(rand(5), checks=false), SO3Algebra)

    # Test construction, with vectors that are not 3-dimensional with checks turned on.
    @test_throws ArgumentError SO3Algebra(rand(2))
    @test_throws ArgumentError SO3Algebra(rand(4))
    @test_throws ArgumentError SO3Algebra(rand(5))

    # Test construction with 3-dimensional vectors.
    @test isa(SO3Algebra([1; 0; 0]), SO3Algebra)
    @test isa(SO3Algebra(rand(3)), SO3Algebra)
end

@testset "SO3Algebra Zero" begin
    test_alg = SO3Algebra(rand(3))

    @test OH + OH ≈ OH
    @test OH + test_alg ≈ test_alg
    @test test_alg + OH ≈ test_alg
end

@testset "SO3Algebra Summation" begin
    @test ALG_X + ALG_X ≈ SO3Algebra([2; 0; 0])
    @test ALG_X + ALG_Y ≈ SO3Algebra([1; 1; 0])
    @test ALG_X + ALG_Z ≈ SO3Algebra([1; 0; 1])
    @test ALG_Y + ALG_X ≈ SO3Algebra([1; 1; 0])
    @test ALG_Y + ALG_Y ≈ SO3Algebra([0; 2; 0])
    @test ALG_Y + ALG_Z ≈ SO3Algebra([0; 1; 1])

    @test ALG_X + ALG_Y + ALG_Z ≈ SO3Algebra([1; 1; 1])
end

@testset "SO3Algebra Negation" begin
    @test -OH ≈ OH
    @test OH + (-OH) ≈ OH
    @test -OH + OH ≈ OH

    @test -ALG_X + ALG_X ≈ OH
    @test -ALG_Y + ALG_Y ≈ OH
    @test -ALG_Z + ALG_Z ≈ OH
    @test ALG_X + (-ALG_X) ≈ OH
    @test ALG_Y + (-ALG_Y) ≈ OH
    @test ALG_Z + (-ALG_Z) ≈ OH
end

@testset "SO3Algebra Subtraction" begin
    @test ALG_X - ALG_X ≈ OH
    @test ALG_X - ALG_Y ≈ SO3Algebra([1; -1; 0])
    @test ALG_X - ALG_Z ≈ SO3Algebra([1; 0; -1])
    @test ALG_Y - ALG_X ≈ SO3Algebra([-1; 1; 0])
    @test ALG_Y - ALG_Y ≈ OH
    @test ALG_Y - ALG_Z ≈ SO3Algebra([0; 1; -1])

    @test ALG_X - ALG_Y - ALG_Z ≈ SO3Algebra([1; -1; -1])
end

@testset "SO3Algebra Logarithmic Map" begin
    @test exp(OH, RotationMatrix) ≈ one(RotationMatrix)
    @test exp(SO3Algebra([0.5 * π; 0; 0]), RotationMatrix) ≈
        RotationMatrix([1 0 0; 0 0 -1; 0 1 0])
    @test exp(SO3Algebra([0; 0.5 * π; 0]), RotationMatrix) ≈
        RotationMatrix([0 0 1; 0 1 0; -1 0 0])
    @test exp(SO3Algebra([0; 0; 0.5 * π]), RotationMatrix) ≈
        RotationMatrix([0 -1 0; 1 0 0; 0 0 1])
end

@testset "SO3Algebra Conversion" begin
    @test convert(Array, OH) ≈ OH.value
end

@testset "SO3Algebra Correctness Checks" begin
    @test check(OH)[1]
    @test check(convert(Array, OH))[1]
    @test !check(SO3Algebra(rand(2), checks=false))[1]
    @test !check(rand(2))[1]
end

@testset "SO3Algebra Approximation Check" begin
    @test OH ≈ OH 
    @test OH ≉ SO3Algebra(rand(3))
end
