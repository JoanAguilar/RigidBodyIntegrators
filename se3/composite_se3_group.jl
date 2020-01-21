"""
    CompositeSE3Group{TypeTransArray<:Real, TypeRot<:AbstractSO3Group} <: AbstractSE3Group

Concrete type that stores an SE3 Group element as a translational part, as a
3-dimensional vector (of type `Array{TypeTransArray, 1}`); and a rotational part (of type
`TypeRot`).

See the `AbstractSE3Gproup` type definition for more information.
"""
struct CompositeSE3Group{TypeTransArray<:Real, TypeRot<:AbstractSO3Group} <:
        AbstractSE3Group
    trans::Array{TypeTransArray, 1}
    rot::TypeRot


    """
        CompositeSE3Group(trans::Array{TypeTransArray, 1}, rot::TypeRot; checks::Bool=true)
	where {TypeTransArray<:Real, TypeRot<:AbstractSO3Group}

    Construct a `CompositeSE3Group` instance from a 3-dimensional vector and a SO3 group
    element.

    `trans` is an 3-dimensional vector; `rot` must be an instance of a subtype of
    `AbstractSO3Group`; if `checks` is set to `true`, correctness checks are performed
    during construction.
    """
    function CompositeSE3Group(
            trans::Array{TypeTransArray, 1}, rot::TypeRot; checks::Bool=true) where
            {TypeTransArray<:Real, TypeRot<:AbstractSO3Group}

        if checks
            # Perform correctness checks.
            if size(trans) ≠ (3,)
                checks_pass = false
                msg = "Expected size (3,) for a translational vector, got $(size(trans)). " *
                    "Argument is not a valid SE3 representation."
            else
                checks_pass, msg = check(rot)
            end
	    if !checks_pass
                throw(ArgumentError(msg))
            end
        end
	new{TypeTransArray, TypeRot}(trans, rot)
    end
end


"""
    Base.:one(in::CompositeSO3Group)

Identity element.
"""
function Base.:one(in::T) where T<:CompositeSE3Group
    return T([0; 0; 0], one(rot(in)), checks=false)
end


# FIXME: It ought to be possible to create a type-based `one` function for
# `CompositeSE3Group`. The trickery is in determining the type of the rotational part of
# the transformation. The implementation below will usually fail with a `MethodError`
# with calls such as `one(CompositeSE3Group{Int64, RotationMatrix})`, due to the element
# type of the `RotationMatrix` type not being specified.
# """
#     Base.:one(type::Type{CompositeSE3Group})
# 
# Identity element of type `type`.
# """
# function Base.:one(type::Type{CompositeSE3Group{TypeArrTrans, TypeRot}}) where
#         {TypeArrTrans<:Real, TypeRot<:AbstractSO3Group}
#     return type([0; 0; 0], one(TypeRot), checks=false)
# end


"""
    Base.:*(left::CompositeSE3Group, right::CompositeSE3Group)

SE3 elements multiplication.

The output `out::CompositeSE3Group`, satisfies the equation `out = left * right`.
"""
function Base.:*(left::CompositeSE3Group, right::CompositeSE3Group)
    trans_new = rot(left) * trans(right) + trans(left)
    rot_new = rot(left) * rot(right)
    return CompositeSE3Group(trans_new, rot_new, checks=false)
end


"""
    Base.:*(gr::CompositeSE3Group, vec::Array{<:Real, 1})

Multiplication by a vector.

The output `out::Array{<:Real, 1}` satisfies `size(out) == (3,)` and is the result of
applying `gr` to `vec`.
"""
function Base.:*(gr::CompositeSE3Group, vec::Array{<:Real, 1})
    return rot(gr) * vec + trans(gr)
end


"""
    Base.:*(gr::CompositeSE3Group, mat::Array{<:Real, 2})

Multiplication by a matrix.

The output `out::Array{<:Real, 2}` satisfies `size(out) == size(mat)` and is the result
of applying `gr` to each of the columns of `mat`.
"""
function Base.:*(gr::CompositeSE3Group, mat::Array{<:Real, 2})
    return rot(gr) * mat + trans(gr) * ones(eltype(trans(gr)), (1, size(mat, 2)))
end


"""
    Base.inv(gr::CompositeSE3Group)

Compute the inverse.

The output `out::CompositeSE3Group` satisfies `out * gr ≈ one(gr)` and
`gr * out ≈ one(gr)` to a certain degree of accuracy.
"""
function Base.inv(gr::CompositeSE3Group)
    return CompositeSE3Group( inv(rot(gr)) * (-trans(gr)), inv(rot(gr)), checks=false)
end


"""
    log(in::CompositeSE3Group)

Logarithmic map.

Convert a `CompositeSE3Group` instance to a `SE3Algebra` instance.
"""
function log(in::CompositeSE3Group)
    return SE3Algebra(trans(in), log(rot(in)), checks=false)
end


"""
    Base.:convert(type::Type{CompositeSE3Group}, gr::CompositeSE3Group; checks::Bool=true)

Convert `gr` to the type `type`. Correctness checks are performed if `checks==true`.
"""
function Base.:convert(type::Type{CompositeSE3Group}, gr::CompositeSE3Group;
                       checks::Bool=true)
    return type(trans(gr), rot(gr), checks=checks)
end


"""
    check(gr::CompositeSE3Group)
 
Check if `gr` is a correct representation of SE3.

Return `true` and an empty string if `gr` is a correct representation of SE3, returns
`false` and an information message, otherwise.
"""
function check(gr::CompositeSE3Group)
    if size(trans(gr)) ≠ (3,)
        msg = "Expected size (3,) for a translational vector, got $(size(trans(gr))). " *
	    "Argument is not a valid SE3 representation."
	return false, msg
    end
    return check(rot(gr))
end

"""
    trans(gr::CompositeSE3Group)

Return the translational part of `gr` as a 3-dimensional vector.
"""
function trans(gr::CompositeSE3Group)
    return gr.trans
end


"""
    rot(gr::CompositeSE3Group)

Return the rotational part of `gr`.
"""
function rot(gr::CompositeSE3Group)
    return gr.rot
end


"""
    rot(type::Type{AbstractSO3Group}, gr::CompositeSE3Group)

Return the rotational part of `gr` as type `type`.
"""
function rot(type::Type{AbstractSO3Group}, gr::CompositeSE3Group)
    return convert(type, gr.rot) 
end


"""
    Base.isapprox(left::RotationMatrix, right::RotationMatrix, args...)

Compare two `RotationMatrix` instances.

See `Base`'s documentation for more information about `args`.
"""
function Base.:isapprox(left::CompositeSE3Group, right::CompositeSE3Group, args...)
    return isapprox(trans(left), trans(right), args...) & isapprox(rot(left), rot(right), args...)
end
