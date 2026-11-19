abstract type ZetaFrontend end

function polynomial end
function cohomology_quotient end
function weight end
function check_input end
function reduction_prescale end
function reduction_ruv_function end
function reduction_directions end
function frobenius_series end
function frobenius_term_pole_order end

function reduction_prescale(P::CohomologyQuotient, m::Integer)
    m >= 1 || throw(ArgumentError("the pole order must be positive"))
    scale = QQ(1)
    for j = 1:(m-1)
        scale *= inv(reduction_scalar(P, j))
    end
    return scale
end

"""
    reduction_ruv_function(fe, S, f, pseudoInverseMat, cache, params)

Return the function that computes the R_uv matrices for each reduction direction.
"""
function reduction_ruv_function(fe::ZetaFrontend, S, f, pseudoInverseMat, cache, params)
    P = cohomology_quotient(fe)
    return V -> computeRuvS(V, S, f, P, pseudoInverseMat, cache, params)
end

"""
    reduction_directions(fe, S)

Return a predicate `(u, k, v) -> Bool` for membership in the direction collection
V_S(u, k) of `fe`.
"""
function reduction_directions(fe::ZetaFrontend, S)
    P = cohomology_quotient(fe)
    d = P.d
    return function (u, k, v)
        length(u) == P.n + 1 || return false
        length(v) == length(u) || return false
        k >= 1 || return false
        sum(v) == d || return false
        all(x -> x >= 0, v) || return false
        all(i -> u[i] - k * v[i] >= 0, eachindex(u)) || return false
        all(i -> 0 <= i <= P.n && (u[i+1] == 0 || v[i+1] != 0), S)
    end
end
