function _quadric_gram_matrix(f, K, n)
    iszero(f) && throw(ArgumentError("expected a nonzero quadratic form"))
    A = Oscar.zero_matrix(K, n, n)
    diagonal = true
    half = inv(K(2))
    for (c, e) in Oscar.coefficients_and_exponents(f)
        sum(e) == 2 || throw(ArgumentError("expected a homogeneous quadratic form"))
        i, j = findfirst(!iszero, e), findlast(!iszero, e)
        diagonal &= i == j
        A[i, j] = A[j, i] = i == j ? c : c * half
    end
    return A, diagonal
end

function _quadric_pencil_data(f, g)
    parent(f) === parent(g) || throw(ArgumentError("expected forms in the same ring"))
    R = parent(f)
    K = Oscar.base_ring(R)
    K isa Oscar.FinField || throw(ArgumentError("expected a finite field"))
    isodd(Oscar.characteristic(K)) || throw(ArgumentError("expected odd characteristic"))
    n = Oscar.nvars(R)
    n >= 3 || throw(ArgumentError("expected at least three variables"))
    A, diagonal_f = _quadric_gram_matrix(f, K, n)
    B, diagonal_g = _quadric_gram_matrix(g, K, n)
    return A, B, K, n, diagonal_f && diagonal_g
end

# At the root [-b_i:a_i] the nonzero block has entries a_i*b_j - b_i*a_j.
# Their product is the cofactor. A zero cross product detects a repeated
# root or a common zero diagonal entry (identically singular determinant).
function _quadric_diagonal_signs(A, B, n)
    a = [A[i, i] for i = 1:n]
    b = [B[i, i] for i = 1:n]
    negative = 0
    for i = 1:n
        cofactor = one(a[i]) * (-1)^((n - 1) ÷ 2)
        for j = 1:n
            i == j && continue
            cofactor *= a[i] * b[j] - b[i] * a[j]
        end
        iszero(cofactor) && return nothing
        negative += !Oscar.is_square(cofactor)
    end
    return negative
end

# For invertible A, det(tA+B) = det(A) charpoly(-A^-1 B)(t).
# If B or B-A is invertible instead, reverse the fixed-length characteristic
# polynomial to obtain det(C+vA), then substitute v=t or v=t+1.
function _quadric_pencil_det(A, B, S, t, d = Oscar.det(A))
    if !iszero(d)
        return d * Oscar.charpoly(S, -inv(A) * B)
    end
    n = Oscar.nrows(A)
    e = Oscar.det(B)
    if !iszero(e)
        return e * reverse(Oscar.charpoly(S, -inv(B) * A), n + 1)
    end
    C = B - A
    e = Oscar.det(C)
    if !iszero(e)
        h = reverse(Oscar.charpoly(S, -inv(C) * A), n + 1)
        return e * h(t + 1)
    end
    M = Oscar.matrix(S, n, n, [t * A[i, j] + B[i, j] for i = 1:n for j = 1:n])
    return Oscar.det(M)
end

# The pivot columns give coordinate vectors complementary to the radical.
function _quadric_pencil_correction(M, q, n)
    k, echelon = Oscar.rref(M)
    isodd(k) && return zero(q)
    k == 0 && return q^(n - 2)
    pivots = Int[findfirst(j -> !iszero(echelon[i, j]), 1:Oscar.ncols(M)) for i = 1:k]
    disc = (-1)^(k ÷ 2) * Oscar.det(M[pivots, pivots])
    return (Oscar.is_square(disc) ? 1 : -1) * q^(n - 2 - k ÷ 2)
end

"""
    quadric_pencil_point_count(f, g)

Count the common projective zeros of two nonzero quadratic forms over an odd
finite field. Counts singular intersections as well. For an odd number of variables, regular members have zero
correction. Singular pencils and even-dimensional Gram matrices are also supported.
"""
function quadric_pencil_point_count(f::MPolyRingElem, g::MPolyRingElem)
    A, B, K, n, diagonal = _quadric_pencil_data(f, g)
    q = BigInt(Oscar.order(K))
    count = sum(q^i for i = 0:(n-3))
    if diagonal && isodd(n)
        negative = _quadric_diagonal_signs(A, B, n)
        if !isnothing(negative)
            return count + (n - 2 * negative) * q^((n - 3) ÷ 2)
        end
    end
    count += _quadric_pencil_correction(A, q, n)
    if q > max(64, n^2)
        leading = Oscar.det(A)
        if iszero(leading) && iszero(Oscar.det(B))
            # The pivot columns span a coordinate complement to ker(A) ∩ ker(B).
            # All ranks and regular-block discriminants survive this restriction.
            k, E = Oscar.rref(vcat(A, B))
            if k < n
                keep = [findfirst(j -> !iszero(E[i, j]), 1:n) for i = 1:k]
                A, B = A[keep, keep], B[keep, keep]
                leading = Oscar.det(A)
            end
        end
        S, t = Oscar.polynomial_ring(K, "t")
        D = _quadric_pencil_det(A, B, S, t, leading)
        if !iszero(D)
            m = Oscar.ncols(A)
            if isodd(m)
                for root in Oscar.roots(D)
                    count += _quadric_pencil_correction(root * A + B, q, n)
                end
            else
                signed_det = (-1)^(m ÷ 2) * D
                balance = zero(q)
                for u in K
                    value = signed_det(u)
                    if iszero(value)
                        count += _quadric_pencil_correction(u * A + B, q, n)
                    else
                        balance += Oscar.is_square(value) ? 1 : -1
                    end
                end
                count += balance * q^(n - 2 - m ÷ 2)
            end
            return count
        end
    end
    for t in K
        count += _quadric_pencil_correction(t * A + B, q, n)
    end
    return count
end

"""
    quadric_pencil_zeta(f, g)

Return the full zeta function of the smooth projective intersection of two
quadrics in an odd number of variables over an odd finite field. A squarefree
binary pencil determinant of full homogeneous degree certifies smoothness.
Raises an `ArgumentError` when the input falls outside this domain.
"""
function quadric_pencil_zeta(f::MPolyRingElem, g::MPolyRingElem)
    A, B, K, n, diagonal = _quadric_pencil_data(f, g)
    isodd(n) || throw(ArgumentError("the zeta formula requires an odd number of variables"))
    R, T = Oscar.polynomial_ring(ZZ, "T")
    q = ZZ(Oscar.order(K))
    if diagonal
        negative = _quadric_diagonal_signs(A, B, n)
        isnothing(negative) && throw(
            ArgumentError("the binary pencil determinant must have only simple roots"),
        )
        scale = q^((n - 3) ÷ 2)
        ambient = prod(1 - q^i * T for i = 0:(n-3))
        return one(R) // (ambient * (1 - scale*T)^(n-negative) * (1 + scale*T)^negative)
    end
    S, t = Oscar.polynomial_ring(K, "t")
    D = _quadric_pencil_det(A, B, S, t)
    if iszero(D) || Oscar.degree(D) < n - 1 || !isone(gcd(D, Oscar.derivative(D)))
        throw(ArgumentError("the binary pencil determinant must have only simple roots"))
    end

    denom = prod(1 - q^i * T for i = 0:(n-3))
    scale = q^((n - 3) ÷ 2)
    sign = (-1)^((n - 1) ÷ 2)
    # At infinity the coefficient of t^(n-1) is the cofactor of A.
    # Corank one guarantees a nonzero principal cofactor in each residue field.
    minors = Dict{Int,typeof(D)}()
    function root_sign(residue_norm)
        for i = 1:n
            minor = get!(minors, i) do
                indices = [j for j = 1:n if j != i]
                sign * _quadric_pencil_det(A[indices, indices], B[indices, indices], S, t)
            end
            value = residue_norm(minor)
            if !iszero(value)
                return Oscar.is_square(value) ? 1 : -1
            end
        end
        error("no nonsingular principal minor at a simple determinant root")
    end
    for (factor, _) in Oscar.factor(D)
        h = Oscar.divexact(factor, Oscar.leading_coefficient(factor))
        epsilon = root_sign(m -> Oscar.resultant(h, m))
        denom *= 1 - epsilon * (scale * T)^Oscar.degree(h)
    end
    if Oscar.degree(D) == n - 1
        denom *= 1 - root_sign(m -> Oscar.coeff(m, n - 1)) * scale * T
    end
    return one(R) // denom
end
