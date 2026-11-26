@testset "pairs of quadrics" begin
    count_pair = DeRham.quadric_pencil_point_count
    zeta_pair = DeRham.quadric_pencil_zeta
    P, T = polynomial_ring(ZZ, "T")

    # On the conic [s^2:s*u:u^2] the second quadric cuts out
    # s^4 + s*u^3 + 2*u^4, irreducible over F_3. No point lies at infinity.
    K = GF(3)
    R, (x, y, z) = polynomial_ring(K, ["x", "y", "z"])
    f, g = x * z - y^2, x^2 + y * z + 2 * z^2
    @test count_pair(f, g) == 0
    @test zeta_pair(f, g) == one(P) // (1 - T^4)

    L, a = Oscar.finite_field(3, 2, "a")
    _, (u, v, w) = polynomial_ring(L, ["u", "v", "w"])
    f2, g2 = u * w - v^2, u^2 + v * w + 2 * w^2
    @test count_pair(f2, g2 + a * f2) == 0
    @test zeta_pair(f2, g2 + a * f2) == one(P) // (1 - T^2)^2

    # P^4 over F_7. Direct enumeration gives 64 common zeros. All five
    # singular pencil members are rational, hence three + signs and two -.
    _, X = polynomial_ring(GF(7), ["x$i" for i = 1:5])
    A = sum(X .^ 2)
    B = sum((i - 1) * X[i]^2 for i = 1:5)
    expected = one(P) // ((1 - T) * (1 - 7T)^4 * (1 + 7T)^2 * (1 - 49T))
    @test count_pair(A, B) == 64
    @test zeta_pair(A, B) == expected
    @test zeta_pair(B, A) == expected  # simple root in the infinity chart
    # Swap coordinates by a shear and change the generators of the pencil.
    Y = [X[1] + X[2]; X[2:end]]
    C, D = sum(Y .^ 2), sum((i - 1) * Y[i]^2 for i = 1:5)
    @test zeta_pair(D, C + D) == expected
    @test count_pair(C, D) == 64

    # In seven variables the discriminant is t^7-t and all seven signs
    # are +1 (each signed cofactor equals the negative derivative, +1).
    _, W = polynomial_ring(GF(7), ["w$i" for i = 1:7])
    ambient = prod(1 - ZZ(7)^i * T for i = 0:4)
    @test zeta_pair(sum(W .^ 2), sum((i - 1) * W[i]^2 for i = 1:7)) ==
          one(P) // (ambient * (1 - 49T)^7)

    # The sole common point is [0:0:1]. Interleaved cross terms exercise
    # rings whose display order differs from the internal term order.
    m = 2x^2 + x*z + y*z
    h = 2x^2 + x*y + 2x*z + y^2 + y*z
    @test count_pair(m, h) == 1

    @test count_pair(x^2, x^2) == 4
    @test count_pair(x * y, x * z) == 5  # a line and an isolated point
    @test_throws ArgumentError zeta_pair(x^2, x^2)
    @test_throws ArgumentError zeta_pair(x^2 + y^2 + z^2, x^2)
    @test_throws ArgumentError zeta_pair(x^2, x^2 + y^2 + z^2)
    @test_throws ArgumentError count_pair(zero(R), f)
    @test_throws ArgumentError count_pair(f + x, g)
    @test_throws ArgumentError count_pair(f, f2)
    _, U = polynomial_ring(GF(2), ["u", "v", "w"])
    @test_throws ArgumentError count_pair(U[1]^2, U[2]^2)
    _, V = polynomial_ring(GF(3), ["v$i" for i = 1:4])
    @test count_pair(V[1]^2, V[2]^2) == 4  # underlying projective line
    @test_throws ArgumentError zeta_pair(V[1]^2, V[2]^2)
end

@testset "rational singular members of odd pencils" begin
    K = GF(101)
    _, (x, y, z) = polynomial_ring(K, ["x", "y", "z"])
    # Repeated determinant root, and a simple root at infinity when swapped.
    # z=0 and (x+y)^2+y^2=0 give two points since -1 is a square.
    f, g = (x + y)^2 + y^2, z^2
    @test DeRham.quadric_pencil_point_count(f, g) == 2
    @test DeRham.quadric_pencil_point_count(g, f) == 2
    @test DeRham.quadric_pencil_point_count(x*y, x^2) == 102

    # The irreducible quartic of the F_3 conic example splits over F_81.
    L, a = Oscar.finite_field(3, 4, "a")
    _, (u, v, w) = polynomial_ring(L, ["u", "v", "w"])
    h, j = u*w - v^2, u^2 + v*w + 2*w^2
    @test DeRham.quadric_pencil_point_count(j + a*h, h) == 4

    # Check diagonal coefficients outside the prime field and a root at infinity.
    diagonal, second = u^2 + v^2 + w^2, a*v^2 + w^2
    mixed, other = (u+v)^2 + v^2 + w^2, a*v^2 + w^2
    @test DeRham.quadric_pencil_zeta(second, diagonal) ==
          DeRham.quadric_pencil_zeta(other, mixed)
    @test DeRham.quadric_pencil_point_count(second, diagonal) ==
          DeRham.quadric_pencil_point_count(other, mixed)
end

@testset "even-rank pencil corrections" begin
    K = GF(101)
    _, X = polynomial_ring(K, ["x$i" for i = 1:4])
    # Direct projective enumeration gives 120 points on the base.
    @test DeRham.quadric_pencil_point_count(sum(X .^ 2), sum((i-1)*X[i]^2 for i = 1:4)) ==
          120
    _, Y = polynomial_ring(K, ["y$i" for i = 1:5])
    U = [Y[i] + Y[5] for i = 1:4]
    # The sheared cone has one vertex and 101 affine points over each base point.
    @test DeRham.quadric_pencil_point_count(sum(U .^ 2), sum((i-1)*U[i]^2 for i = 1:4)) ==
          12121
    _, (x, y, z) = polynomial_ring(K, ["x", "y", "z"])
    # Every member is singular, but the two radicals have zero intersection.
    @test DeRham.quadric_pencil_point_count(x*z, y*z) == 103
end
