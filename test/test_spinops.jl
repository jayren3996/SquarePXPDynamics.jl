@testset "spin operators" begin
    X = pauli_x()
    Y = pauli_y()
    Z = pauli_z()
    I2 = identity2()

    @test X * X ≈ I2
    @test Y * Y ≈ I2
    @test Z * Z ≈ I2
    @test X * Y ≈ im * Z
    @test projector_up() ≈ ComplexF64[1 0; 0 0]
    @test projector_down() ≈ ComplexF64[0 0; 0 1]
    @test projector_up() + projector_down() ≈ I2
    @test projector_up() * projector_down() ≈ zeros(ComplexF64, 2, 2)
    @test Z * projector_up() ≈ projector_up()
    @test Z * projector_down() ≈ -projector_down()
    @test kron_all([X, Z]) ≈ kron(X, Z)
    @test_throws ArgumentError kron_all(Matrix{ComplexF64}[])
    @test embed_one_site(X, 1, 2) ≈ kron(X, I2)
    @test embed_one_site(X, 2, 2) ≈ kron(I2, X)
    @test size(embed_one_site(X, 1, 3)) == (8, 8)
    @test_throws ArgumentError embed_one_site(X, 0, 2)
    @test_throws ArgumentError embed_one_site(ones(ComplexF64, 3, 3), 1, 2)
end
