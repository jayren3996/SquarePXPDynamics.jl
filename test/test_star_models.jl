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

@testset "PXP ITensor star gate convention" begin
    sites = [Index(2, tag) for tag in ("center", "right", "up", "left", "down")]
    model = PXPStarModel(false)
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
    model = PXPStarModel(true)
    protocol = StaticModel(model)
    @test model_at(protocol, 0.0, 1) === model
    @test model_at(protocol, 1.0, 17) === model
end
