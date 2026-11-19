@testset "preconfigured R_uv rebinding" begin
    coefficient_ring, = residue_ring(ZZ, 5)
    matspace = matrix_space(coefficient_ring, 1, 1)
    first_matrix = matspace()
    second_matrix = matspace()
    first_matrix[1, 1] = 1
    second_matrix[1, 1] = 2

    pep = DeRham.PregenLazyPEP{typeof(first_matrix)}(_ -> [first_matrix])
    point = [1]
    original_matrix = only(pep[point])

    DeRham.bind_ruv_function!(pep, _ -> [second_matrix])

    @test only(pep[point]) === original_matrix
    @test original_matrix[1, 1] == 2
end
