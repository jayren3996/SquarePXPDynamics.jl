@testset "spin operators" begin
    X = pauli_x()
    Y = pauli_y()
    Z = pauli_z()
    I2 = identity2()

    @test X * X ≈ I2
    @test Y * Y ≈ I2
    @test Z * Z ≈ I2
    @test X * Y ≈ im * Z
    @test projector_up() + projector_down() ≈ I2
    @test projector_up() * projector_down() ≈ zeros(ComplexF64, 2, 2)
end
