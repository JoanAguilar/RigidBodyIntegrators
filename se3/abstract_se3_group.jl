include("../abstract_lie.jl")

"""
AbstractSE3Group

Abstract type representing SE3 elements.

# Implementation
Subtypes of of this type (`AbstractSE3Group`) are implementations of representations of SE3
group elements (referred as `<Implementation>`). The expectation is for them to implement:
 - `Base.:one`: Identity element. Satisfies `I * elem == elem`, `elem * I == elem`, for
   any `elem::<Implementation>`.
 - `Base.:*(left::<Implementation>, right::<Implementation>)`: Multiplication between
   elements. The result is an element representing the application of `right` and then
   `left`.
 - `Base.:*(elem::<Implementation>, vec::Array{<:Real, 1})`: Multiplication between an
   element and a 3-dimensional vector. The result is a vector obtained by applying `elem`
   to `vec`.  
 - `Base.:*(elem::<Implementation>, mat::Array{<:Real, 2})`: Multiplication between an
   element and a 3×n matrix. The result is a 3×n matrix obtained by applying `elem` to
   each of the columns of `mat`.
 - `Base.:inv(elem::<Implementation>)`: Element inverse. `inv` satisfies that
   `elem * inv(elem)` or `inv(elem) * elem` returns the identity element `one`.
 - `log(elem::<Implementation>)`: Logarithmic map. The result is an instance of
   `SE3Algebra` equivalent to `elem`.
 - `Base.:convert(impl::Type{AbstractSO3Group}, elem::<Implementation>)`: Convert the
   element `elem` to type `impl`.
 - `check(elem::<Implementation>)`: Check for correctness. Returns `true` and an empty
   string if `elem` is a correct representation of SO3. Returns `false` and an
   information message, otherwise.
 - `trans(elem::<Implementation>)`: Return the translational part of `elem` as a
   3-dimensional vector.
 - `rot(elem::<Implementation>; type::Type{AbstractSO3Group})`: Return the translational
   part of `elem`. If `type` is specified, the output satisfies `elem::type`.
"""
abstract type AbstractSE3Group <: AbstractLieGroup end
