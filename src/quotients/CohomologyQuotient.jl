"""
    abstract type CohomologyQuotient

A Griffiths–Dwork quotient model of cohomology indexed by pole order.
"""
abstract type CohomologyQuotient end

"""
    HypersurfaceQuotient(n, d)

The cohomology quotient for a degree-`d` hypersurface of dimension `n`.
"""
struct HypersurfaceQuotient <: CohomologyQuotient
    n::Int
    d::Int
end

"""
    pole_orders(P::HypersurfaceQuotient)

Return the pole orders that carry basis elements for the hypersurface `P`.
"""
pole_orders(P::HypersurfaceQuotient) = 1:P.n

"""
    terminal_pole_order(P::HypersurfaceQuotient)

Return the pole order where reduction of the hypersurface `P` stops.
"""
terminal_pole_order(P::HypersurfaceQuotient) = P.n

"""
    numerator_degree(P::HypersurfaceQuotient, m)

Return the degree of the numerator at pole order `m` for `P`.
"""
numerator_degree(P::HypersurfaceQuotient, m::Integer) = P.d * m - (P.n + 1)

"""
    pole_exponent(P::HypersurfaceQuotient, m)

Return the exponent of the hypersurface polynomial at pole order `m`.
"""
pole_exponent(::HypersurfaceQuotient, m::Integer) = m

"""
    reduction_scalar(P::HypersurfaceQuotient, m)

Return the exact Griffiths–Dwork scalar for the step from pole order `m + 1` to `m`.
"""
reduction_scalar(::HypersurfaceQuotient, m::Integer) = QQ(1, m)

"""
    DoubleCoverQuotient(n, d)

The cohomology quotient for the double cover of `P^n` branched along a
degree-`d` hypersurface, with `d` even.
"""
struct DoubleCoverQuotient <: CohomologyQuotient
    n::Int
    d::Int

    function DoubleCoverQuotient(n::Integer, d::Integer)
        iseven(d) || throw(ArgumentError("the branch degree must be even"))
        new(Int(n), Int(d))
    end
end

"""
    pole_orders(P::DoubleCoverQuotient)

Return the pole orders that carry basis elements for the double cover `P`.
"""
pole_orders(P::DoubleCoverQuotient) = 1:(P.n+1)

"""
    terminal_pole_order(P::DoubleCoverQuotient)

Return the pole order where reduction of the double cover `P` stops.
"""
terminal_pole_order(P::DoubleCoverQuotient) = P.n + 1

"""
    numerator_degree(P::DoubleCoverQuotient, m)

Return the degree of the numerator at pole order `m` for `P`.
"""
numerator_degree(P::DoubleCoverQuotient, m::Integer) = P.d * m - (P.n + 1 + P.d ÷ 2)

"""
    pole_exponent(P::DoubleCoverQuotient, m)

Return the exponent in the double cover pole `y^(2m-1)`.
"""
pole_exponent(::DoubleCoverQuotient, m::Integer) = 2m - 1

"""
    reduction_scalar(P::DoubleCoverQuotient, m)

Return the exact Griffiths–Dwork scalar for the step from pole order `m + 1`
to `m`.
"""
reduction_scalar(::DoubleCoverQuotient, m::Integer) = QQ(2, 2m - 1)
