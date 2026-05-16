using ITensors
using LinearAlgebra

function _assert_sane_simple_gauge_diag(diag; near_zero::Bool)
    @test diag isa SimpleGaugeDiagnostic
    @test diag.bond isa BondKey
    @test isfinite(diag.deviation)
    @test diag.deviation >= 0
    @test isfinite(diag.frobenius_norm)
    @test diag.frobenius_norm > 0
    @test isfinite(diag.diagonal_min)
    @test isfinite(diag.diagonal_max)
    @test isfinite(diag.diagonal_sum)
    @test diag.diagonal_min >= -1e-12
    @test diag.diagonal_max >= diag.diagonal_min
    @test diag.diagonal_sum >= diag.diagonal_max - 1e-12
    if near_zero
        @test diag.deviation ≈ 0 atol = 1e-12
    end
end

function _diagonal_product_square_ipeps(cell)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for c in cell.reps
        p = physical_index(psi, c)
        left = link_index(psi, c, :left)
        right = link_index(psi, c, :right)
        up = link_index(psi, c, :up)
        down = link_index(psi, c, :down)
        T = ITensor(ComplexF64, p, left, right, up, down)
        for a = 1:2
            T[p=>2, left=>a, right=>a, up=>a, down=>a] = 1.0 + 0.0im
        end
        psi.tensors[c] = T
        set_link_weight!(psi, c, :right, [0.8, 0.6])
        set_link_weight!(psi, c, :up, [0.7, 0.5])
    end
    return psi
end

function _seeded_offdiagonal_square_ipeps(cell)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for c in cell.reps
        p = physical_index(psi, c)
        left = link_index(psi, c, :left)
        right = link_index(psi, c, :right)
        up = link_index(psi, c, :up)
        down = link_index(psi, c, :down)
        T = ITensor(ComplexF64, p, left, right, up, down)
        for pv = 1:dim(p), lv = 1:dim(left), rv = 1:dim(right), uv = 1:dim(up), dv = 1:dim(down)
            re = 0.01 * (13c.x + 17c.y + 5pv + 3lv + 7rv + 11uv + dv)
            im = 0.001 * (2pv + lv - rv + uv - dv)
            T[p=>pv, left=>lv, right=>rv, up=>uv, down=>dv] = complex(re, im)
        end
        psi.tensors[c] = T
        set_link_weight!(psi, c, :right, [0.9, 0.4])
        set_link_weight!(psi, c, :up, [0.6, 0.3])
    end
    return psi
end

function _copy_with_gauge(psi, gauge::Symbol)
    return SquareIPEPSState(
        psi.unitcell,
        copy(psi.tensors),
        copy(psi.physical_indices),
        copy(psi.link_indices),
        Dict(key => copy(values) for (key, values) in psi.link_weights),
        psi.maxdim,
        gauge,
        Ref(state_version(psi)),
        Ref(log_norm(psi)),
    )
end

@testset "simple gauge diagnostics D=1 product state" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(2, 2)

    right = gauge_diagnostic_simple(psi, c, :right)
    left_alias = gauge_diagnostic_simple(psi, neighbor(cell, c, :right), :left)
    @test right.bond == bondkey(cell, c, :right)
    @test left_alias.bond == right.bond
    @test left_alias.deviation ≈ right.deviation atol = 1e-14
    @test left_alias.frobenius_norm ≈ right.frobenius_norm atol = 1e-14
    @test gauge_deviation_simple(psi, c, :right) ≈ right.deviation atol = 1e-14
    _assert_sane_simple_gauge_diag(right; near_zero = true)
    _assert_sane_simple_gauge_diag(left_alias; near_zero = true)

    deviations = all_gauge_deviations_simple(psi)
    @test deviations isa Dict{BondKey,Float64}
    @test length(deviations) == length(psi.link_weights)
    @test Set(keys(deviations)) == Set(keys(psi.link_weights))
    @test all(isfinite, values(deviations))
    @test all(value -> isapprox(value, 0; atol = 1e-12), values(deviations))
end

@testset "simple gauge diagnostics D=2 diagonal product state" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = _diagonal_product_square_ipeps(cell)
    c = SquareCoord(1, 3)

    for dir in (:right, :up)
        diag = gauge_diagnostic_simple(psi, c, dir)
        @test diag.bond == bondkey(cell, c, dir)
        @test gauge_deviation_simple(psi, c, dir) ≈ diag.deviation atol = 1e-14
        _assert_sane_simple_gauge_diag(diag; near_zero = true)
    end
end

@testset "simple gauge diagnostics detect off-diagonal D=2 fixture" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = _seeded_offdiagonal_square_ipeps(cell)
    deviations = all_gauge_deviations_simple(psi)

    @test length(deviations) == length(psi.link_weights)
    @test all(isfinite, values(deviations))
    @test any(value -> value > 1e-10, values(deviations))

    diag = gauge_diagnostic_simple(psi, SquareCoord(2, 1), :right)
    _assert_sane_simple_gauge_diag(diag; near_zero = false)
    @test diag.deviation > 1e-10
end

@testset "simple gauge diagnostics validation" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    c = SquareCoord(1, 1)

    @test_throws ArgumentError gauge_diagnostic_simple(psi, c, :diagonal)
    @test_throws ArgumentError gauge_deviation_simple(psi, c, :diagonal)

    not_simple = _copy_with_gauge(psi, :not_simple)
    @test not_simple.gauge === :not_simple
    @test psi.gauge === :simple
    @test not_simple.unitcell === psi.unitcell
    @test not_simple.maxdim == psi.maxdim
    @test state_version(not_simple) == state_version(psi)
    @test log_norm(not_simple) == log_norm(psi)
    @test not_simple !== psi
    @test not_simple.link_weights !== psi.link_weights
    @test all(key -> not_simple.link_weights[key] == psi.link_weights[key], keys(psi.link_weights))

    @test_throws ArgumentError gauge_diagnostic_simple(not_simple, c, :right)
    @test_throws ArgumentError gauge_deviation_simple(not_simple, c, :right)
    @test_throws ArgumentError all_gauge_deviations_simple(not_simple)
end
