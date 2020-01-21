using LinearAlgebra: norm

"""
    SO3Algebra{T<:Real}

Concrete type that stores an SO3 Algebra element as 3-element vector.
"""
struct SO3Algebra{T<:Real} <: AbstractLieAlgebra
    value::Array{T, 1}

    """
        SO3Algebra(value::Array{T, 1}; checks::Bool=true) where T<:Real

    Construct a `SO3Algebra` instance from a vector (`value`).

    If `checks` is set to `true`, correctness checks are performed during
    construction.
    """
    function SO3Algebra(value::Array{T, 1}; checks::Bool=true) where T<:Real
        if checks
            # Perform correctness checks.
            checks_pass, msg = check_as(value, SO3Algebra)
            if !checks_pass
                throw(ArgumentError(msg))
            end
        end
        new{T}(value)
    end
end


"""
    Base.:zero(in::SO3Algebra)

Zero element.
"""
function Base.:zero(in::T) where T<:SO3Algebra
    return T([0; 0; 0], checks=false)
end

"""
    Base.:zero(type::Type{SO3Algebra})

Zero element of type `type`.
"""
function Base.:zero(type::Type{SO3Algebra})
    return type([0; 0; 0], checks=false)
end


"""
    Base.:+(left::SO3Algebra, right::SO3algebra)

Summation.
"""
function Base.:+(left::SO3Algebra, right::SO3Algebra)
    return SO3Algebra(left.value + right.value, checks=false)
end


"""
    Base.:-(in::SO3Algebra)

Negation.

The output `out` satisfies `in + out = zero`.
"""
function Base.:-(in::SO3Algebra)
    return SO3Algebra(-in.value, checks=false)
end


"""
    Base.:-(left::SO3Algebra, right::SO3Algebra)

Subtraction.
"""
function Base.:-(left::SO3Algebra, right::SO3Algebra)
    return SO3Algebra(left.value - right.value, checks=false)
end


"""
    exp(in::SO3Algebra, type::Type{RotationMatrix})

Exponential map.

The output `out` will satisfy `out::type`.
"""
function exp(in::SO3Algebra, type::Type{RotationMatrix})
    # NOTE: It might be possible to have a "lazy" version of this function
    # (`exp(::SO3Algebra)`), where the output type is unspecified at runtime. This "lazy"
    # version should allow for things like `exp(a) * b`, where the result of `exp(a)` is
    # computed and matched with the type of `b` during the multiplication; or
    # `exp(zero(SO3Algebra)) ≈ one(type)`, where the result of `exp` is computed and
    # matched with the type `type` during comparison.
    in_arr = convert(Array, in)
    θ = norm(in_arr)
    if θ > 0
        ax = in_arr / θ
	k = [0 -ax[3] ax[2];
             ax[3] 0 -ax[1];
             -ax[2] ax[1] 0]
        mat = I + sin(θ) * k + (1 - cos(θ)) * k^2
        return type(mat, checks=true)
    else
        return one(type)
    end
end


"""
    Base.:convert(type::Type{Array}, alg::SO3Algebra; copy::Bool=false)

Convert `alg` to `type`. If `copy == true`, a copy of the array data will be returned.
"""
function Base.:convert(type::Type{Array}, alg::SO3Algebra; copy::Bool=false)
    if copy
        return type(copy(alg.value))
    else
        return type(alg.value)
    end
end


"""
    check(in::SO3Algebra)

Check if `in` is a correct representation of an SO3 algebra element.

Return `true` and an empty string if `in` is a correct representation, returns `false`
and an information message, otherwise.
"""
function check(in::SO3Algebra)
    return check_as(in.value, typeof(in))
end


"""
    check_as(value::Array{<:Real, 1}, type::Type{SO3Algebra})

Check if `value` is a 3-element vector. See `check(in::SO3Algebra)` for more information
on return types.
"""
function check_as(value::Array{<:Real, 1}, type::Type{<:SO3Algebra})
    if size(value) ≠ (3,)
        msg = "Expected size (3,) got $size(value). Argument can't be used as a " *
        "representation of an SO3 algebra element."
        return false, msg
    else
        return true, ""
    end
end


"""
    Base.isapprox(left::SO3Algebra, right::SO3Algebra, args...)

Compare two `SO3Algebra` instances.

See `Base`'s documentation for more information about `args`.
"""
function Base.isapprox(left::SO3Algebra, right::SO3Algebra, args...)
    return isapprox(left.value, right.value, args...)
end
