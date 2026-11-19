function _frontend_direction_chunk_length(u, v, terminal)
    mins = copy(u)
    m = terminal + 4
    k = 0
    while true
        next = mins .- v
        any(next .< 0) && break
        m - k == terminal && break
        mins = next
        k += 1
    end
    return k
end

@testset "frontend reduction directions" begin
    Rcurve, (x, z) = polynomial_ring(GF(5), ["x", "z"])
    curve = DeRham.ProjectiveHypersurface(x^3 + z^3)
    Rsurface, (a, b, c) = polynomial_ring(GF(5), ["a", "b", "c"])
    surface = DeRham.ProjectiveHypersurface(a^3 + b^3 + c^3)
    Rfourfold, variables = polynomial_ring(GF(5), ["x$i" for i = 0:4])
    fourfold = DeRham.ProjectiveHypersurface(sum(variable^3 for variable in variables))
    Rdouble_curve, (s, t) = polynomial_ring(GF(5), ["s", "t"])
    double_curve = DeRham.DoubleCover(s^4 + t^4)
    Rdouble_surface, (u, v, w) = polynomial_ring(GF(5), ["u", "v", "w"])
    double_surface = DeRham.DoubleCover(u^6 + v^6 + w^6)

    depthfirst_cases = (
        (curve, [0, 1], [3, 3]),
        (surface, [0, 1, 2], [3, 3, 3]),
        (fourfold, [0], [3, 0, 0, 0, 0]),
        (double_curve, [0, 1], [4, 4]),
        (double_surface, [0, 1, 2], [6, 6, 6]),
    )
    for (fe, S, u) in depthfirst_cases
        P = DeRham.cohomology_quotient(fe)
        contains_direction = DeRham.reduction_directions(fe, S)
        V = DeRham.rev_chooseV(copy(u), P.d, S)
        k = _frontend_direction_chunk_length(u, V, DeRham.terminal_pole_order(P))
        @test k > 0
        @test contains_direction(u, k, V)
    end

    varbyvar_cases = (
        (curve, [1], [1, 4]),
        (surface, [2], [1, 1, 4]),
        (fourfold, [4], [1, 1, 1, 1, 3]),
    )
    for (fe, S, u) in varbyvar_cases
        P = DeRham.cohomology_quotient(fe)
        contains_direction = DeRham.reduction_directions(fe, S)
        V = DeRham.varbyvar_chooseV(copy(u), P.d)
        k = _frontend_direction_chunk_length(u, V, DeRham.terminal_pole_order(P))
        @test k > 0
        @test contains_direction(u, k, V)
    end

    contains_first_variable = DeRham.reduction_directions(curve, [0])
    @test !contains_first_variable([3, 3], 1, [0, 3])
end
