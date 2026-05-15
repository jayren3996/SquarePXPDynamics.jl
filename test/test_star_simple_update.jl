using ITensors
using LinearAlgebra

function _star_coords_test(cell, center)
    c = wrap(cell, center)
    return (
        center = c,
        right = neighbor(cell, c, :right),
        up = neighbor(cell, c, :up),
        left = neighbor(cell, c, :left),
        down = neighbor(cell, c, :down),
    )
end

function _dense_square_star_index_test(values)
    idx = 1
    for (site, value) in enumerate(values)
        idx += (value - 1) * 2^(SQUARE_STAR_SITES - site)
    end
    return idx
end

function _site_amplitudes_d1(psi, c)
    p = physical_index(psi, c)
    T = psi.tensors[wrap(psi.unitcell, c)]
    others = filter(!=(p), inds(T))
    return [T[p => value, (i => 1 for i in others)...] for value in 1:2]
end

function _local_density_d1_testhelper(psi, c)
    amp = _site_amplitudes_d1(psi, c)
    normsq = sum(abs2, amp)
    normsq > 0 || throw(ArgumentError("zero D=1 site norm"))
    return real(abs2(amp[1]) / normsq)
end

function _dense_star_state_d1_testhelper(psi, center)
    coords = _star_coords_test(psi.unitcell, center)
    sites = (coords.center, coords.right, coords.up, coords.left, coords.down)
    amplitudes = [_site_amplitudes_d1(psi, c) for c in sites]
    state = zeros(ComplexF64, 2^SQUARE_STAR_SITES)
    for values in Iterators.product((1:2 for _ in 1:SQUARE_STAR_SITES)...)
        idx = _dense_square_star_index_test(values)
        state[idx] = prod(amplitudes[i][values[i]] for i in 1:SQUARE_STAR_SITES)
    end
    return state
end

function _assert_finite_star_update_state_test(psi, infos = StarUpdateInfo[])
    for T in values(psi.tensors)
        @test isfinite(norm(T))
    end
    for lambda in values(psi.link_weights)
        @test all(isfinite, lambda)
        @test norm(lambda) ≈ 1 atol = 1e-12
    end
    for info in infos
        @test isfinite(info.max_truncerr)
        @test all(isfinite, values(info.truncerrs))
        @test all(isfinite, values(info.norm_factors))
        @test all(isfinite, values(info.min_lambda))
    end
end

function _normalized_dense_star_state_d1_testhelper(psi, center)
    state = _dense_star_state_d1_testhelper(psi, center)
    return state ./ norm(state)
end

function _assert_dense_states_equivalent_test(actual, expected; atol = 1e-10)
    pivot = argmax(abs.(expected))
    phase = expected[pivot] / actual[pivot]
    @test norm(phase .* actual - expected) < atol
end

function _star_site_set_test(cell, centers)
    sites = Set{SquareCoord}()
    for center in centers
        coords = _star_coords_test(cell, center)
        for site in (coords.center, coords.right, coords.up, coords.left, coords.down)
            push!(sites, site)
        end
    end
    return sites
end

function _bond_endpoint_set_test(cell, bond)
    c = bond.site
    other = neighbor(cell, c, bond.dir)
    return Set((c, other))
end

@testset "star validation" begin
    psi_bad = product_square_ipeps(PeriodicSquareUnitCell(2, 2); state = :down, maxdim = 1)
    @test_throws ArgumentError project_star!(psi_bad, SquareCoord(1, 1), 0.1)

    psi = product_square_ipeps(PeriodicSquareUnitCell(10, 10); state = :down, maxdim = 1)
    info = project_star!(psi, SquareCoord(5, 5), 0.0; maxdim = 1)
    @test info isa StarUpdateInfo
end

@testset "invalid split order" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(10, 10); state = :down, maxdim = 1)
    @test_throws ArgumentError project_star!(
        psi,
        SquareCoord(5, 5),
        0.1;
        split_order = (:right, :up, :left, :left),
    )
end

@testset "zero-step star update preserves D=1 checkerboard diagnostics" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
    density_before = density_simple(psi)
    blockade_before = blockade_violation_simple(psi)
    weights_before = deepcopy(psi.link_weights)

    info = project_star!(psi, SquareCoord(5, 5), 0.0; maxdim = 1)

    @test info.max_truncerr ≈ 0 atol = 1e-12
    @test density_simple(psi) ≈ density_before atol = 1e-12
    @test blockade_violation_simple(psi) ≈ blockade_before atol = 1e-12
    for (bond, lambda) in psi.link_weights
        @test norm(lambda) ≈ 1 atol = 1e-12
        @test length(lambda) == 1
        @test lambda[1] ≈ weights_before[bond][1] atol = 1e-12
    end
end

@testset "all-down flippable center update" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    dt = 0.05

    info = project_star!(psi, c, dt; evolution = :real, projected = true, maxdim = 1)

    ncenter = _local_density_d1_testhelper(psi, c)
    @test ncenter ≈ sin(dt)^2 atol = 1e-10 rtol = 1e-8
    @test info.max_truncerr ≈ 0 atol = 1e-12
end

@testset "D=1 star update matches dense local reference" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    dt = 0.07

    before = _dense_star_state_d1_testhelper(psi, c)
    expected = projected_square_pxp_gate(dt; evolution = :real) * before
    project_star!(psi, c, dt; evolution = :real, projected = true, maxdim = 1)
    actual = _dense_star_state_d1_testhelper(psi, c)

    expected ./= norm(expected)
    actual ./= norm(actual)
    _assert_dense_states_equivalent_test(actual, expected)
end

@testset "star update diagnostics finite" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    info = project_star!(psi, SquareCoord(5, 5), 0.03; maxdim = 1)

    dirs = Set([:right, :up, :left, :down])
    @test isfinite(info.max_truncerr)
    @test Set(keys(info.truncerrs)) == dirs
    @test Set(keys(info.keptdims)) == dirs
    @test Set(keys(info.norm_factors)) == dirs
    @test Set(keys(info.min_lambda)) == dirs
    @test all(isfinite, values(info.truncerrs))
    @test all(>=(0), values(info.truncerrs))
    @test all(k -> k <= 1, values(info.keptdims))
    @test all(isfinite, values(info.norm_factors))
    @test all(isfinite, values(info.min_lambda))
end

@testset "repeated D=1 all-down star updates stay finite and normalized" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    infos = StarUpdateInfo[]

    for _ in 1:10
        push!(infos, project_star!(psi, c, 0.01; evolution = :real, projected = true, maxdim = 1))
        _assert_finite_star_update_state_test(psi, infos)
    end
end

@testset "same-color disjoint star updates leave far-away link weights unchanged" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    centers = (SquareCoord(5, 5), SquareCoord(10, 5))
    @test square_star_color(centers[1]) == square_star_color(centers[2])
    @test stars_are_disjoint_mod_unitcell(cell, collect(centers))

    weights_before = deepcopy(psi.link_weights)
    affected_sites = _star_site_set_test(cell, centers)

    info1 = project_star!(psi, centers[1], 0.02; evolution = :real, projected = true, maxdim = 1)
    info2 = project_star!(psi, centers[2], 0.02; evolution = :real, projected = true, maxdim = 1)

    _assert_finite_star_update_state_test(psi, [info1, info2])
    for (bond, lambda_before) in weights_before
        isempty(intersect(_bond_endpoint_set_test(cell, bond), affected_sites)) || continue
        @test psi.link_weights[bond] == lambda_before
    end
end

@testset "D=1 non-default split order matches dense star reference" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    dt = 0.07

    psi_default = product_square_ipeps(cell; state = :down, maxdim = 1)
    psi_reordered = product_square_ipeps(cell; state = :down, maxdim = 1)
    before = _dense_star_state_d1_testhelper(psi_default, center)
    expected = projected_square_pxp_gate(dt; evolution = :real) * before
    expected ./= norm(expected)

    project_star!(psi_default, center, dt; evolution = :real, projected = true, maxdim = 1)
    info = project_star!(
        psi_reordered,
        center,
        dt;
        evolution = :real,
        projected = true,
        maxdim = 1,
        split_order = (:up, :right, :down, :left),
    )

    actual_default = _normalized_dense_star_state_d1_testhelper(psi_default, center)
    actual_reordered = _normalized_dense_star_state_d1_testhelper(psi_reordered, center)
    _assert_dense_states_equivalent_test(actual_reordered, expected)
    _assert_dense_states_equivalent_test(actual_reordered, actual_default)
    @test _local_density_d1_testhelper(psi_reordered, center) ≈ sin(dt)^2 atol = 1e-10 rtol = 1e-8
    _assert_finite_star_update_state_test(psi_reordered, [info])
end

@testset "invalid update after valid update does not partially corrupt state" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(10, 10); state = :down, maxdim = 1)
    project_star!(psi, SquareCoord(5, 5), 0.03; evolution = :real, projected = true, maxdim = 1)
    weights_before = deepcopy(psi.link_weights)
    density_before = density_simple(psi)
    blockade_before = blockade_violation_simple(psi)

    @test_throws ArgumentError project_star!(
        psi,
        SquareCoord(5, 5),
        0.1;
        split_order = (:right, :right, :up, :down),
    )

    @test psi.link_weights == weights_before
    @test density_simple(psi) ≈ density_before atol = 1e-12
    @test blockade_violation_simple(psi) ≈ blockade_before atol = 1e-12
end

@testset "D=2 star update smoke test" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)

    info = project_star!(psi, SquareCoord(5, 5), 0.03; maxdim = 2)

    @test isfinite(info.max_truncerr)
    @test all(isfinite, values(info.truncerrs))
    @test all(k -> k <= 2, values(info.keptdims))
    for lambda in values(psi.link_weights)
        @test isfinite(norm(lambda))
        @test norm(lambda) ≈ 1 atol = 1e-12
    end
end
