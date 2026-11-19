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

@testset "frontend Frobenius term conversion" begin
    R, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    fe = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    P = DeRham.cohomology_quotient(fe)
    precision_ring, = residue_ring(ZZ, 25)
    Rlift, (x_lift, y_lift, z_lift) =
        polynomial_ring(precision_ring, ["x", "y", "z"])
    S = [0, 1, 2]
    params = DeRham.ZetaFunctionParams(
        0,
        false,
        :depthfirst,
        :invlex,
        false,
        true,
        false,
        false,
        false,
    )
    cache = DeRham.controlled_reduction_cache(P.n, P.d, S, params, P)
    order = 10
    lifted_polynomial = 2 * x_lift^3 * z_lift^2 + 3 * y_lift^3 * z_lift^2
    keys = DeRham.reduction_input_keys(fe, [[lifted_polynomial, order]], S)
    @test length(keys) == 1
    @test keys[1][1] == order
    @test length(keys[1][2]) == 1
    input = keys[1][2][1]
    @test input.key == ([0, 1, 3], order)

    basis = cache[DeRham._reduction_g_degree(P)]
    g = zeros(UInt, length(basis))
    DeRham.reduction_input_vector!(g, fe, input, S, cache, params)
    expected = zeros(UInt, length(g))
    expected[findfirst(==([4, 0, 0]), basis)] = 2
    expected[findfirst(==([1, 3, 0]), basis)] = 3
    @test g == expected
    @test count(x -> !iszero(x), g) == 2

    hypersurface_output = [
        [x^3 * z^2, 5],
        [y^3 * z^2, 10],
        [z^5, 15],
    ]
    hypersurface_orders = [
        group[1] for group in DeRham.reduction_input_keys(fe, hypersurface_output, S)
    ]
    @test sort(hypersurface_orders) == sort(collect(15:-5:5))

    Rdouble, (s, t) = polynomial_ring(GF(5), ["s", "t"])
    double_cover = DeRham.DoubleCover(s^4 + t^4)
    double_cover_output = [[s^4, 3], [t^4, 8], [s^4, 13]]
    double_cover_orders = [
        group[1] for group in DeRham.reduction_input_keys(double_cover, double_cover_output, [1])
    ]
    @test double_cover_orders == [3, 8, 13]
    @test 3 in double_cover_orders

    terminal_S = [0]
    terminal_cache = DeRham.controlled_reduction_cache(P.n, P.d, terminal_S, params, P)
    basis_exponents = terminal_cache[DeRham._reduction_g_degree(P)]
    exponent_index = findfirst(==([4, 0, 0]), basis_exponents)
    terminal_vector = zeros(UInt, length(basis_exponents))
    terminal_vector[exponent_index] = 1
    terminal = DeRham.reduction_terminal_polynomial(
        fe,
        [([1, 0, 0], terminal_vector)],
        parent(lifted_polynomial),
        terminal_S,
        params,
    )
    @test terminal == [[x_lift^4, DeRham.terminal_pole_order(P)]]
end
