include("abstract_so3_group.jl")
include("so3_algebra.jl")

using LinearAlgebra: det


"""
    RotationMatrix{T<:Real} <: AbstractSO3Group

Concrete type that stores an SO3 Group element as a 3×3 rotation matrix.

See the `AbstractSO3Gproup` type definition for more information.
"""
struct RotationMatrix{T<:Real} <: AbstractSO3Group
    # Rotation matrix.
    value::Array{T, 2}


    """
        RotationMatrix(value::Array{T, 2}; checks::Bool=true) where {T<:Real}

    Construct a `RotationMatrix` instance from an array.

    `value` is an array with the 3×3 rotation matrix; if `checks` is set to `true`,
    correctness checks are performed during construction.
    """
    function RotationMatrix(value::Array{T, 2}; checks::Bool=true) where {T<:Real}
        if checks
            # Perform correctness checks.
            checks_pass, msg = check(value) 
	    if !checks_pass
                throw(ArgumentError(msg))
            end
        end
	new{T}(value)
    end
end


"""
    one::RotationMatrix{Int64}

Identity element.
"""
one = RotationMatrix([1 0 0; 0 1 0; 0 0 1], checks=false)


"""
    Base.:*(left::RotationMatrix, right::RotationMatrix)

Rotation matrix multiplication.

The output `out::RotationMatrix`, satisfies the equation `out = left * right`.
"""
function Base.:*(left::RotationMatrix, right::RotationMatrix)
    return RotationMatrix(left.value * right.value, checks=false)
end


"""
    Base.:*(mat::RotationMatrix, vec::Array{<:Real, 1})

Multiplication by a vector.

The output `out::Array{<:Real, 1}` satisfies `size(out) == (3,)` and the equation
`out = mat * vec`.
"""
function Base.:*(mat::RotationMatrix, vec::Array{<:Real, 1})
    return mat.value * vec
end


"""
    Base.:*(rotmat::RotationMatrix, mat::Array{<:Real, 2})

Multiplication by a matrix.

The output `out::Array{<:Real, 2}` satisfies `size(out) == size(mat)` and the equation
`out = rotmat * mat`.
"""
function Base.:*(rotmat::RotationMatrix, mat::Array{<:Real, 2})
    return rotmat.value * mat
end


"""
    Base.inv(mat::RotationMatrix)

Compute the inverse.

The output `out::RotationMatrix` satisfies `out * mat ≈ I` and `mat * out ≈ I` to a
certain degree of accuracy.
"""
function Base.inv(mat::RotationMatrix)
    return RotationMatrix(convert(Array, mat.value'), checks=false)
end


"""
    log(mat::RotationMatrix)

Logarithmic map.

Convert a `RotationMatrix` instance to a `SO3Algebra` instance.
"""
function log(mat::RotationMatrix)
    ax = mat.value
    θ = acos((mat[1, 1] + mat[2, 2] + mat[3, 3] - 1) / 2)
    if θ > 0
        ax_den = 2 * sin(θ)
        ax_num_vec = [mat[3, 2] - mat[2, 3],
                      mat[1, 3] - mat[3, 1],
                      mat[2, 1] - mat[1, 2]]
        return SO3Algebra(θ * ax_num_vec / ax_den, checks=false)
    else
        return SO3Algebra(zeros(eltype(mat), 3), checks=false)
    end
end


"""
    convert(type, mat::RotationMatrix; checks::Bool=true)

Convert a `mat` to the type `type`. Correctness checks are performed if `checks==true`.
Note that if the output `out` satisfies `out<:RotationMatrix` the array information is
copied.
"""
function convert(type, mat::RotationMatrix; checks::Bool=true)
    if type<:RotationMatrix
        return T(copy(mat.value), checks=checks)
    else
        return T(mat.value, checks=checks)
    end
end


"""
    check(mat::RotationMatrix)

Check if `mat` is a correct representation of SO3.

Return `true` and an empty string if `mat` is a correct representation of SO3, returns
`false` and an information message, otherwise.
"""
function check(mat::RotationMatrix)
    return check(mat.value)
end


"""
    check(value::Array{<:Real, 2})

Check if `value` is a valid rotation matrix. See `check(mat::RotationMatrix)` for more
information on the return types.
"""
function check(value::Array{<:Real, 2})
    if size(value) ≠ (3, 3)
        # Check that the given matrix is 3×3.
        msg = "Expected size (3, 3), got $size(value). Argument is not a rotation matrix."
	return false, msg
    end
    if value * value' ≉ LinearAlgebra.I
        # Check that the transpose is the inverse.
        msg = "Transpose is not an inverse. Argument is not a rotation matrix."
	return false, msg
    end
    if det(value) ≉ 1
        # Check that the determinant is 1.
        msg = "Expected determinant 1, got $det(value). Argument is not a rotation matrix."
	return false, msg
    end
    return true, ""
end


"""
    asarray(mat::RotationMatrix; copy::Bool=false)

Return an array with the rotation matrix. If `copy == true`, a copy of the array data will
be returned.
"""
function asarray(mat::RotationMatrix; copy::Bool=false)
    if copy
        return copy(mat.value)
    else
        return mat.value
    end
end


"""
    Base.isapprox(left::RotationMatrix, right::RotationMatrix, args...)

Compare two `RotationMatrix` instances.

See `Base`'s documentation for more information about `args`.
"""
function Base.isapprox(left::RotationMatrix, right::RotationMatrix, args...)
    return isapprox(left.value, right.value, args...)
end
