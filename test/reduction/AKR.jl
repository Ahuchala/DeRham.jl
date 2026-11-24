@testset "AKR algorithm identifiers" begin
    @test DeRham.akr === :akr
    @test DeRham.akr_naive === :akr_naive
    @test DeRham.akr_levelbylevel === :akr_levelbylevel

    R, (x, y, z) = polynomial_ring(GF(5), ["x", "y", "z"])
    hypersurface = DeRham.ProjectiveHypersurface(x^3 + y^3 + z^3)
    for algorithm in (DeRham.akr, DeRham.akr_naive, DeRham.akr_levelbylevel)
        @test DeRham._supports_reduction_algorithm(hypersurface, algorithm)
    end

end

@testset "AKR plane curve full-precision equality" begin
    curve = DeRham.fermat_hypersurface(3, 3, 2)
    _, naive =
        DeRham.zeta_coefficients(curve; givefrobmat = true, algorithm = DeRham.akr_naive)
    @test naive == BigInt[2, 0, 1]
    @test DeRham.zeta_coefficients(curve; algorithm = DeRham.akr) == naive
    @test DeRham.zeta_coefficients(curve; algorithm = DeRham.akr_levelbylevel) == naive
    @test DeRham.zeta_coefficients(curve) == naive
end

@testset "AKR quartic K3 limited-precision matrices" begin
    _, (x, y, z, w) = polynomial_ring(GF(2), ["x", "y", "z", "w"])
    k3 = x^3*y + y^3*z + z^3*w + w^3*x + x*y*z*w
    default_matrix = DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = :default)
    naive_matrix =
        DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = DeRham.akr_naive)
    levelbylevel_matrix =
        DeRham.frobenius_matrix_with_precision(k3, 1; algorithm = DeRham.akr_levelbylevel)
    @test naive_matrix == default_matrix
    @test levelbylevel_matrix == default_matrix
end
