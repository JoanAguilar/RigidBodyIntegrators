"""
    SO3Algebra{T<:Real}

Concrete type that stores an SO3 Algebra element as 3-element vector.
"""
struct SO3Algebra{T<:Real}
    value::Array{T, 1}

    """
        SO3Algebra(value::Array{T, 1}; checks::Bool=true) where {T<:Real}

    Construct a `SO3Algebra` instance from a vector (`value`).

    If `checks` is set to `true`, correctness checks are performed during
    construction.
    """
    function SO3Algebra(value::Array{T, 1}; checks::Bool=true) where {T<:Real}
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
    Base.:zero(in::SO3Algebra)

Zero element.
"""
function Base.:zero(in::SO3Algebra)
    return typeof(in)([0; 0; 0], checks=false)
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
    exp(in::SO3Algebra, type)

Exponential map.

The output `out` will satisfy `out::type`.
"""
function exp(in::SO3Algebra, type)
    if type<:RotationMatrix
        θ = norm(in.value)
        if θ > 0
            ax = in.value / θ
	    k = skew(ax)
	    mat = I + sin(θ) * k + (1 - cos(θ)) * k
	    return T(mat, checks=true)
        else
            return one
        end
    else
        throw(ErrorException("Can't instantiate exponential map to $T, not implemented."))
    end
end


"""
    check(in::SO3Algebra)

Check if `in` is a correct representation of an SO3 algebra element.

Return `true` and an empty string if `in` is a correct representation, returns `false`
and an information message, otherwise.
"""
function check(in::SO3Algebra)
    return check(in.value)
end


"""
    check(value::Array{<:Real, 1})

Check if `value` is a 3-element vector. See `check(in::SO3Algebra)` for more information
on return types.
"""
function check(value::Array{<:Real, 1})
    if size(value) ≠ (3,)
        msg = "Expected size (3,) got $size(value). Argument can't be used as a " *
	"representation of an SO3 algebra element."
	return false, msg
    else
        return true, ""
    end
end


"""
    asarray(in::SO3Algebra; copy::Bool=false)

Return an array with the 3-element vector associated with `in`. If `copy == true`, a copy
of the vector data will be returned.
"""
function asarray(in::SO3Algebra; copy::Bool=false)
    if copy
        return copy(in.value)
    else
        return in.value
    end
end


"""
    Base.isapprox(left::RotationMatrix, right::RotationMatrix, args...)

Compare two `SO3Algebra` instances.

See `Base`'s documentation for more information about `args`.
"""
function Base.isapprox(left::SO3Algebra, right::SO3Algebra, args...)
    return isapprox(left.value, right.value, args...)
end

