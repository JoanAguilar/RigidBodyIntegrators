"""
AbstractSO3Group

Abstract type representing SO3 elements.

# Implementation
Subtypes of of this type (AbstractSO3Group) are implementations of representations of SO3
group elements (referred as `<Implementation>`). The expectation is for them to implement:
 - `I::<Implementation>`: Identity element. Satisfies `I * elem == elem`,
    `elem * I == elem`, for any `elem::<Implementation>`.
 - `Base.:*(left::<Implementation>, right::<Implementation>)`: Multiplication between
   elements. The result is an element representing the application of `right` and then
   `left`.
 - `Base.:*(elem::<Implementation>, vec::Array{<:Real, 1})`: Multiplication between an
   element and a 3-dimensional vector. The result is a vector obtained by applying `elem`
   to `vec`.  
 - `Base.:*(elem::<Implementation>, mat::Array{<:Real, 2})`: Multiplication between an
   element and a 3×n matrix. The result is a 3×n matrix obtained by applying `elem` to
   each of the columns of `mat`.
 - `inv(elem::<Implementation>)`: Element inverse. `inv` satisfies `elem * inv(elem)` or
   `inv(elem) * elem`. Returns the identity element `I`.
 - `log(elem::<Implementation>)`: Logarithmic map. The result is an instance of
   SO3Algebra equivalent to `elem`.
 - `convert(impl, elem::<Implementation>)`: Convert the element `elem` to type `impl`.
   `impl` must be an existing SO3 Group implementation.
 - `check(elem::<Implementation>)`: Check for correctness. Returns `true` and an empty
   string if `elem` is a correct representation of SO3. Returns `false` and an
   information message, otherwise.
"""
abstract type AbstractSO3Group end
