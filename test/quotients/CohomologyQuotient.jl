@testset "cohomology quotients" begin
    for (n, d) in ((1, 3), (1, 4), (2, 4), (3, 3))
        quotient = DeRham.HypersurfaceQuotient(n, d)
        @test quotient isa DeRham.CohomologyQuotient
        @test DeRham.pole_orders(quotient) == 1:n
        @test DeRham.terminal_pole_order(quotient) == n
        @test DeRham.numerator_degree(quotient, n) == d * n - (n + 1)
        for m = 1:n
            @test DeRham.numerator_degree(quotient, m) == d * m - (n + 1)
            @test DeRham.pole_exponent(quotient, m) == m
            @test DeRham.reduction_scalar(quotient, m) == Oscar.QQ(1, m)
        end
    end

    hypersurface = DeRham.HypersurfaceQuotient(1, 3)
    double_cover = DeRham.DoubleCoverQuotient(1, 4)
    @test [DeRham.reduction_prescale(hypersurface, m) for m = 1:6] ==
          [factorial(ZZ(m - 1)) for m = 1:6]
    @test [DeRham.reduction_prescale(double_cover, m) for m = 1:6] ==
          [QQ(1), QQ(1, 2), QQ(3, 4), QQ(15, 8), QQ(105, 16), QQ(945, 32)]

    @test_throws ArgumentError DeRham.DoubleCoverQuotient(2, 3)

    Rcurve, (x, y) = polynomial_ring(GF(11), ["x", "y"])
    params = DeRham.default_params()
    for genus = 1:4
        d = 2 * genus + 2
        f = x^d + y^d
        P = DeRham.DoubleCoverQuotient(1, d)
        @test P isa DeRham.CohomologyQuotient
        @test DeRham.pole_orders(P) == 1:2
        @test DeRham.terminal_pole_order(P) == 2
        for m in DeRham.pole_orders(P)
            @test DeRham.numerator_degree(P, m) == d * m - (2 + d ÷ 2)
            @test DeRham.pole_exponent(P, m) == 2m - 1
            @test DeRham.reduction_scalar(P, m) == Oscar.QQ(2, 2m - 1)
        end
        cache = DeRham.controlled_reduction_cache(1, d, [0, 1], params, P)
        basis = DeRham.compute_monomial_bases(f, P, params, cache)
        @test length.(basis) == [genus, genus]
        @test sum(length, basis) == 2 * genus
    end

    Rsurface, (a, b, c) = polynomial_ring(GF(7), ["a", "b", "c"])
    sextic = a^6 + b^6 + c^6
    surface_quotient = DeRham.DoubleCoverQuotient(2, 6)
    surface_cache =
        DeRham.controlled_reduction_cache(2, 6, [0, 1, 2], params, surface_quotient)
    surface_basis =
        DeRham.compute_monomial_bases(sextic, surface_quotient, params, surface_cache)
    @test length.(surface_basis) == [1, 19, 1]
    @test sum(length, surface_basis) == 21

    quartic = a^4 + b^4 + c^4
    quartic_quotient = DeRham.DoubleCoverQuotient(2, 4)
    quartic_cache =
        DeRham.controlled_reduction_cache(2, 4, [0, 1, 2], params, quartic_quotient)
    quartic_basis =
        DeRham.compute_monomial_bases(quartic, quartic_quotient, params, quartic_cache)
    @test isempty(quartic_basis[1])
    @test isempty(quartic_basis[end])
end
