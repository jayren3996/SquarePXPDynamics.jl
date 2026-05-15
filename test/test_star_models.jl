using ITensors
using LinearAlgebra

function _star_models_dense_index_test(values)
    idx = 1
    for (site, value) in enumerate(values)
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

function _star_models_gate_entry_test(G, sites, out_values, in_values)
    out = prime.(sites)
    return G[
        (out[i]=>out_values[i] for i in eachindex(sites))...,
        (sites[i]=>in_values[i] for i in eachindex(sites))...,
    ]
end

function _star_models_dense_from_itensor_gate_test(G, sites)
    dense = zeros(ComplexF64, 2^SQUARE_STAR_SITES, 2^SQUARE_STAR_SITES)
    for out_values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
        out_idx = _star_models_dense_index_test(out_values)
        for in_values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
            in_idx = _star_models_dense_index_test(in_values)
            dense[out_idx, in_idx] =
                _star_models_gate_entry_test(G, sites, out_values, in_values)
        end
    end
    return dense
end

@testset "star model conventions" begin
    @test star_site_order() == (:center, :right, :up, :left, :down)
    @test tfim_pauli_convention() == (:Z_up_is_plus_one, :X_field)
end

@testset "PXP star model reproduces existing gates" begin
    dt = 0.037
    @test star_hamiltonian(PXPStarModel(false)) ≈ square_pxp_star_hamiltonian()
    @test star_gate(PXPStarModel(false), dt; evolution = :real) ≈
          square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(true), dt; evolution = :real) ≈
          projected_square_pxp_gate(dt; evolution = :real)
    @test star_gate(PXPStarModel(false), dt; evolution = :imaginary) ≈
          square_pxp_gate(dt; evolution = :imaginary)
    @test star_gate(PXPStarModel(true), dt; evolution = :imaginary) ≈
          projected_square_pxp_gate(dt; evolution = :imaginary)
end

@testset "TFIM dense star Hamiltonian convention" begin
    model = TFIMStarModel(2.0, 3.0)
    H = star_hamiltonian(model)
    @test size(H) == (32, 32)
    @test H ≈ H'
    @test tfim_product_basis_energy(model, (:up, :up, :up, :up, :up)) ≈ -4.0
    @test tfim_product_basis_energy(model, (:up, :down, :up, :down, :up)) ≈ 0.0
    @test tfim_product_basis_energy(model, (:down, :up, :up, :up, :up)) ≈ 4.0
    @test H[
        _star_models_dense_index_test((1, 1, 1, 1, 1)),
        _star_models_dense_index_test((1, 1, 1, 1, 1)),
    ] ≈ tfim_product_basis_energy(model, (:up, :up, :up, :up, :up))
    @test H[
        _star_models_dense_index_test((2, 1, 1, 1, 1)),
        _star_models_dense_index_test((2, 1, 1, 1, 1)),
    ] ≈ tfim_product_basis_energy(model, (:down, :up, :up, :up, :up))
end

@testset "TFIM dense gates" begin
    model = TFIMStarModel(1.0, 0.7)
    dt = 0.011
    U = star_gate(model, dt; evolution = :real)
    G = star_gate(model, dt; evolution = :imaginary)
    @test U' * U ≈ I atol = 1e-12 rtol = 1e-12
    @test all(isfinite, real.(G))
    @test norm(G' * G - I) > 1e-6
    @test star_gate(model, dt; evolution = :real) *
          star_gate(model, dt; evolution = :real) ≈
          star_gate(model, 2dt; evolution = :real) atol = 1e-12 rtol = 1e-12
    @test_throws ArgumentError star_gate(model, dt; evolution = :bad)
    @test_throws ArgumentError star_gate(model, Inf; evolution = :real)
    @test_throws ArgumentError star_gate(model, NaN; evolution = :imaginary)
    @test_throws ArgumentError TFIMStarModel(Inf, 1.0)
    @test_throws ArgumentError TFIMStarModel(1.0, NaN)
end

@testset "TFIM ITensor star gate convention" begin
    sites = [Index(2, tag) for tag in ("center", "right", "up", "left", "down")]
    model = TFIMStarModel(1.0, 0.7)
    dt = 0.011
    G = star_gate_itensor(model, sites, dt; evolution = :real)

    @test inds(G) == (prime.(sites)..., sites...)
    @test _star_models_dense_from_itensor_gate_test(G, sites) ≈
          star_gate(model, dt; evolution = :real)
    @test _star_models_dense_from_itensor_gate_test(
        star_gate_itensor(model, dt, sites; evolution = :imaginary),
        sites,
    ) ≈ star_gate(model, dt; evolution = :imaginary)
    @test_throws ArgumentError star_gate_itensor(model, sites[1:4], dt)
    @test_throws ArgumentError star_gate_itensor(
        model,
        [Index(3, "bad"), sites[2:end]...],
        dt,
    )
end

@testset "static model protocol" begin
    model = TFIMStarModel(1.0, 2.0)
    protocol = StaticModel(model)
    @test model_at(protocol, 0.0, 1) === model
    @test model_at(protocol, 1.0, 17) === model
end
