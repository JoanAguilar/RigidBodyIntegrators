using Random: rand, seed!
using Test: @test, @testset

seed!(0)

# FIXME: The definition below should be using the `one` function.
EYE = CompositeSE3Group([0; 0; 0], one(RotationMatrix))
XZ = CompositeSE3Group(
    [1; 0; 0],
    RotationMatrix([0 -1 0;
                    1  0 0;
                    0  0 1]))
YY = CompositeSE3Group(
    [0; 1; 0],
    RotationMatrix([0 0 1;
                    0 1 0;
                   -1 0 0]))
ZX = CompositeSE3Group(
    [0; 0; 1],
    RotationMatrix([1 0  0;
                    0 0 -1;
                    0 1  0]))


@testset "Construction" begin
    # Test construction, with checks turned off.
    @test isa(CompositeSE3Group(rand(2), one(RotationMatrix), checks=false),
              CompositeSE3Group)
    @test isa(CompositeSE3Group(rand(3), one(RotationMatrix), checks=false),
              CompositeSE3Group)
    @test isa(CompositeSE3Group(rand(4), one(RotationMatrix), checks=false),
              CompositeSE3Group)

    # Test construction, with translation vectors that are not 3-dimensional with checks
    # turned on.
    @test_throws ArgumentError CompositeSE3Group(rand(2), one(RotationMatrix))
    @test_throws ArgumentError CompositeSE3Group(rand(4), one(RotationMatrix))
end

@testset "Identity" begin
    @test EYE * EYE ≈ EYE
    @test EYE * XZ ≈ XZ
    @test EYE * YY ≈ YY
    @test EYE * ZX ≈ ZX
    @test XZ * EYE ≈ XZ
    @test YY * EYE ≈ YY
    @test ZX * EYE ≈ ZX
end

@testset "Multiplication" begin
    # CompositeSE3Group-CompositeSE3Group multiplication
    @test XZ * YY ≈ CompositeSE3Group([0; 0; 0], RotationMatrix([0 -1 0; 0 0 1; -1 0 0]))
    @test XZ * ZX ≈ CompositeSE3Group([1; 0; 1], RotationMatrix([0 0 1; 1 0 0; 0 1 0]))
    @test YY * XZ ≈ CompositeSE3Group([0; 1; -1], RotationMatrix([0 0 1; 1 0 0; 0 1 0]))
    @test YY * ZX ≈ CompositeSE3Group([1; 1; 0], RotationMatrix([0 1 0; 0 0 -1; -1 0 0]))
    @test ZX * XZ ≈ CompositeSE3Group([1; 0; 1], RotationMatrix([0 -1 0; 0 0 -1; 1 0 0]))
    @test ZX * YY ≈ CompositeSE3Group([0; 0; 2], RotationMatrix([0 0 1; 1 0 0; 0 1 0]))

    # CompositeSE3Group-vector multiplication
    v = [1; 2; 3]
    @test XZ * v ≈ [-1; 1; 3]
    @test YY * v ≈ [3; 3; -1]
    @test ZX * v ≈ [1; -3; 3]

    # CompositeSE3Group-matrix multiplication
    mat = [1 2 3;
	   2 3 1;
	   3 1 2]
    @test XZ * mat ≈ [-1 -2 0;
                       1  2 3;
                       3  1 2]
    @test YY * mat ≈ [3  1  2;
                      3  4  2;
                     -1 -2 -3]
    @test ZX * mat ≈ [1  2  3;
                     -3 -1 -2;
                      3  4  2]
end

@testset "Inverse" begin
    @test inv(EYE) ≈ EYE
    @test EYE * inv(EYE) ≈ EYE
    @test inv(EYE) * EYE ≈ EYE

    @test inv(XZ) * XZ ≈ EYE
    @test inv(YY) * YY ≈ EYE
    @test inv(ZX) * ZX ≈ EYE
    @test XZ * inv(XZ) ≈ EYE
    @test YY * inv(YY) ≈ EYE
    @test ZX * inv(ZX) ≈ EYE
end

@testset "Logarithmic Map" begin
    # FIXME: This should actually be:
    # @test log(one(CompositeSE3Group)) ≈ zero(SE3Algebra)
    # but it fails.
    @test log(EYE) ≈ zero(SE3Algebra)
    @test log(XZ) ≈ SE3Algebra([1; 0; 0], [0; 0; 0.5 * π])
    @test log(YY) ≈ SE3Algebra([0; 1; 0], [0; 0.5 * π; 0])
    @test log(ZX) ≈ SE3Algebra([0; 0; 1], [0.5 * π; 0; 0])
end


# TODO: Conversion tests (`Base.:convert`) are missing. To be implemented once more types of
# SE3Group representations are implemented.


@testset "Correctness Checks" begin
    @test check(EYE)[1]
    @test !check(CompositeSE3Group(rand(2), RotationMatrix(rand(3, 3), checks=false),
				   checks=false))[1]
    @test !check(CompositeSE3Group(rand(4), RotationMatrix(rand(3, 3), checks=false),
				   checks=false))[1]
end


@testset "Translational Getter" begin
    @test trans(EYE) ≈ [0; 0; 0]
    @test trans(XZ) ≈ [1; 0;0 ]
    @test trans(YY) ≈ [0; 1; 0]
    @test trans(ZX) ≈ [0; 0; 1]
end


@testset "Rotational Getter" begin
    @test rot(EYE) ≈ one(RotationMatrix)
    @test rot(XZ) ≈ RotationMatrix([0 -1 0; 1 0 0; 0 0 1])
    @test rot(YY) ≈ RotationMatrix([0 0 1; 0 1 0; -1 0 0])
    @test rot(ZX) ≈ RotationMatrix([1 0 0; 0 0 -1; 0 1 0])

    # TODO: Write tests for the rotational getter where the the output is converted to the
    # desired type (`rot(type::Type{AbstractSO3Group}, gr::CompositeSE3Group)`).
end


@testset "Approximation Check" begin
    @test EYE ≈ EYE
    @test EYE ≉ CompositeSE3Group(rand(3), RotationMatrix(rand(3, 3), checks=false), checks=false)
end
