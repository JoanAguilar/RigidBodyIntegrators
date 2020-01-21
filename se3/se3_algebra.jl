"""
    SE3Algebra{T<:Real} <: AbstractLieAlgebra

Concrete type that stores an SE3 Algebra element.
"""
struct SE3Algebra{T<:Real} <: AbstractLieAlgebra
    trans::Array{T, 1}
    rot::SO3Algebra

    """
        SE3Algebra(trans::Array{T, 1}, rot::SO3Algebra; checks::Bool=true) where {T<:Real}

    Construct a `SE3Algebra` instance from a 3-dimensional vector and a `SO3Algebra` element.

    `trans` is an 3-dimensional vector; `rot` must be an instance of `SO3Algebra`; if `checks` is
    set to `true`, correctness checks are performed during construction.
    """
    function SE3Algebra(trans::Array{T, 1}, rot::SO3Algebra; checks::Bool=true) where
            {T<:Real}
        if checks
            if size(trans) ≠ (3,)
                checks_pass = false
                msg = "Expected size (3,), got $(size(value))"
            else
                checks_pass, msg = check(rot)
            end
            if !checks_pass
                throw(ArgumentError(msg))
            end
        end
	new{T}(trans, rot)
    end
end


"""
    SE3Algebra(trans::Array{TypeTrans, 1}, rot::Array{TypeRot}; checks::Bool=true) where
    {TypeTrans<:Real, TypeRot<:Real}

Construct a `SE3Algebra` instance from two 3-dimensional vectors (one for the translational part
and one for the rotational part).

If `checks` is set to `true`, correctness checks are performed during construction.
"""
function SE3Algebra(trans::Array{TypeTrans, 1}, rot::Array{TypeRot, 1}; checks::Bool=true) where
        {TypeTrans<:Real, TypeRot<:Real}
    if checks
        if size(trans) ≠ (3,)
	    throw(ArgumentError("Expected size (3,), got $(size(value))"))
	end
    end
    new_rot = SO3Algebra(rot, checks=checks)
    return SE3Algebra(trans, new_rot)
end


"""
    Base.:zero(in::SE3Algebra)

Zero element.
"""
function Base.:zero(in::T) where T<:SE3Algebra
    return T([0; 0; 0], zero(SO3Algebra), checks=false)
end


"""
    Base.:zero(type::Type{SE3Algebra})

Zero element of type `type`.
"""
function Base.:zero(type::Type{SE3Algebra})
    return type([0; 0; 0], zero(SO3Algebra), checks=false)
end


"""
    Base.:+(left::SE3Algebra, right::SE3algebra)

Summation.
"""
function Base.:+(left::SE3Algebra, right::SE3Algebra)
    return SE3Algebra(trans(left) + trans(right), rot(left) + rot(right), checks=false)
end


"""
    Base.:-(in::SE3Algebra)

Negation.

The output `out` satisfies `in + out = zero`.
"""
function Base.:-(in::SE3Algebra)
    return SE3Algebra(-trans(in), -rot(in), checks=false)
end


"""
    Base.:-(left::SE3Algebra, right::SE3Algebra)

Subtraction.
"""
function Base.:-(left::SE3Algebra, right::SE3Algebra)
    return SE3Algebra(trans(left) - trans(right), rot(left) - rot(right), checks=false)
end


"""
    exp(in::SE3Algebra, type::Type{CompositeSE3Group})

Exponential map.

The output `out` will satisfy `out::type`.
"""
function exp(in::SE3Algebra, type::Type{CompositeSE3Group})
    return type(trans(in), exp(rot(in)), checks=false)
end


"""
   Base.:convert(type::Type{Array}, alg::SE3Algebra; copy::Bool=false)

Convert `alg` to `type`. If `copy == true`, a copy of the array data will be returned.
"""
function Base.:convert(type::Type{Array}, alg::SE3Algebra; copy::Bool=false)
    arr = type(hcat(trans(alg), convert(type, rot(alg))))
    if copy
        return copy(arr)
    else
        return arr
    end
end


"""
    check(in::SE3Algebra)

Check if `in` is a correct representation of an SE3 algebra element.

Return `true` and an empty string if `in` is a correct representation, returns `false`
and an information message, otherwise.
"""
function check(in::SE3Algebra)
    if size(trans(in)) ≠ (3,)
        return false, "Expected size (3,), got $(size(value))"
    else
        return check(rot(in))
    end
end


"""
    trans(alg::SE3Algebra)

Return the translational part of `alg` as a 3-dimensional vector.
"""
function trans(alg::SE3Algebra)
    return alg.trans
end


"""
    rot(alg::SE3Algebra)

Return the rotational part of `alg`.
"""
function rot(alg::SE3Algebra)
    return alg.rot
end


"""
    Base.isapprox(left::SE3Algebra, right::SE3Algebra, args...)

Compare two `SE3Algebra` instances.

See `Base`'s documentation for more information about `args`.
"""
function Base.:isapprox(left::SE3Algebra, right::SE3Algebra, args...)
    return isapprox(trans(left), trans(right), args...) &
        isapprox(rot(left), rot(right), args...)
end
