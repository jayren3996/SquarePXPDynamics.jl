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
    return [T[p=>value, (i=>1 for i in others)...] for value = 1:2]
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
    for values in Iterators.product((1:2 for _ = 1:SQUARE_STAR_SITES)...)
        idx = _dense_square_star_index_test(values)
        state[idx] = prod(amplitudes[i][values[i]] for i = 1:SQUARE_STAR_SITES)
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

function _assert_finite_simple_summary_test(summary)
    @test all(
        isfinite,
        (
            summary.density,
            summary.density_even,
            summary.density_odd,
            summary.blockade_violation,
            summary.pxp_energy_density,
            summary.mean_bond_entropy,
            summary.max_bond_entropy,
        ),
    )
end

function _seeded_nontrivial_d2_ipeps_test(cell)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for c in cell.reps
        p = physical_index(psi, c)
        left = link_index(psi, c, :left)
        right = link_index(psi, c, :right)
        up = link_index(psi, c, :up)
        down = link_index(psi, c, :down)
        T = ITensor(ComplexF64, p, left, right, up, down)
        for pv = 1:dim(p), lv = 1:dim(left), rv = 1:dim(right), uv = 1:dim(up), dv = 1:dim(down)
            re = 0.01 * (11c.x + 7c.y + 5pv + 3lv + 2rv + uv + dv)
            T[p=>pv, left=>lv, right=>rv, up=>uv, down=>dv] = complex(re, 0.0)
        end
        psi.tensors[c] = T
        set_link_weight!(psi, c, :right, [0.8, 0.6])
        set_link_weight!(psi, c, :up, [0.6, 0.8])
    end
    return psi
end

function _has_active_second_virtual_sector_test(psi)
    for (c, T) in psi.tensors
        p = physical_index(psi, c)
        left = link_index(psi, c, :left)
        right = link_index(psi, c, :right)
        up = link_index(psi, c, :up)
        down = link_index(psi, c, :down)
        for pv = 1:dim(p), lv = 1:dim(left), rv = 1:dim(right), uv = 1:dim(up), dv = 1:dim(down)
            if (lv == 2 || rv == 2 || uv == 2 || dv == 2) &&
               abs(T[p=>pv, left=>lv, right=>rv, up=>uv, down=>dv]) > 0
                return true
            end
        end
    end
    return false
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
    @test state_version(psi) == 1

    @test_throws ArgumentError project_star!(psi, SquareCoord(5, 5), Inf)
    @test_throws ArgumentError project_star!(psi, SquareCoord(5, 5), NaN)
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

@testset "star update records pre-update touched link minima" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    center = SquareCoord(5, 5)
    right_leaf = neighbor(cell, center, :right)
    up_leaf = neighbor(cell, center, :up)

    set_link_weight!(psi, center, :right, [0.3, 0.9])
    set_link_weight!(psi, right_leaf, :right, [0.2, 0.8])
    set_link_weight!(psi, up_leaf, :up, [0.4, 0.7])

    info = project_star!(psi, center, 0.01; evolution = :real, projected = true, maxdim = 2)

    internal = bondkey(cell, center, :right)
    external_right = bondkey(cell, right_leaf, :right)
    external_up = bondkey(cell, up_leaf, :up)
    @test info.touched_min_lambda[internal] ≈ 0.3
    @test info.touched_min_lambda[external_right] ≈ 0.2
    @test info.touched_min_lambda[external_up] ≈ 0.4
    @test length(info.touched_min_lambda) == 16
    @test all(key -> key isa BondKey, keys(info.touched_min_lambda))
    @test all(isfinite, values(info.touched_min_lambda))
    @test all(>=(0), values(info.touched_min_lambda))
end

@testset "repeated D=1 all-down star updates stay finite and normalized" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    c = SquareCoord(5, 5)
    infos = StarUpdateInfo[]

    for _ = 1:10
        push!(
            infos,
            project_star!(psi, c, 0.01; evolution = :real, projected = true, maxdim = 1),
        )
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

    info1 = project_star!(
        psi,
        centers[1],
        0.02;
        evolution = :real,
        projected = true,
        maxdim = 1,
    )
    info2 = project_star!(
        psi,
        centers[2],
        0.02;
        evolution = :real,
        projected = true,
        maxdim = 1,
    )

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
    @test _local_density_d1_testhelper(psi_reordered, center) ≈ sin(dt)^2 atol = 1e-10 rtol =
        1e-8
    _assert_finite_star_update_state_test(psi_reordered, [info])
end

@testset "invalid update after valid update does not partially corrupt state" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(10, 10); state = :down, maxdim = 1)
    project_star!(
        psi,
        SquareCoord(5, 5),
        0.03;
        evolution = :real,
        projected = true,
        maxdim = 1,
    )
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
    @test psi.maxdim == 2

    info = project_star!(psi, SquareCoord(5, 5), 0.03; maxdim = 3)

    @test isfinite(info.max_truncerr)
    @test all(isfinite, values(info.truncerrs))
    @test all(k -> k <= 3, values(info.keptdims))
    @test psi.maxdim == 2
    for lambda in values(psi.link_weights)
        @test isfinite(norm(lambda))
        @test norm(lambda) ≈ 1 atol = 1e-12
    end
end

@testset "repeated D=2 star updates stay finite and normalized" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    c = SquareCoord(5, 5)
    infos = StarUpdateInfo[]

    for _ = 1:3
        push!(
            infos,
            project_star!(psi, c, 0.01; evolution = :real, projected = true, maxdim = 2),
        )
        _assert_finite_star_update_state_test(psi, infos)
        summary = measure_simple(psi)
        @test all(
            isfinite,
            (
                summary.density,
                summary.blockade_violation,
                summary.pxp_energy_density,
                summary.mean_bond_entropy,
                summary.max_bond_entropy,
            ),
        )
    end
end

@testset "D=2 zero-step and split-order observable regressions" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    Hstar = square_pxp_star_hamiltonian()
    psi_zero = product_square_ipeps(cell; state = :down, maxdim = 2)
    density_before = local_density_simple(psi_zero, center)
    star_before = star_expectation_simple(psi_zero, center, Hstar)
    log_norm_before = log_norm(psi_zero)

    project_star!(psi_zero, center, 0.0; projected = false, maxdim = 2)

    @test local_density_simple(psi_zero, center) ≈ density_before atol = 1e-12
    @test star_expectation_simple(psi_zero, center, Hstar) ≈ star_before atol = 1e-12
    @test log_norm(psi_zero) ≈ log_norm_before atol = 1e-12

    psi_default = product_square_ipeps(cell; state = :down, maxdim = 2)
    psi_reordered = product_square_ipeps(cell; state = :down, maxdim = 2)
    project_star!(psi_default, center, 0.03; evolution = :real, projected = true, maxdim = 2)
    project_star!(
        psi_reordered,
        center,
        0.03;
        evolution = :real,
        projected = true,
        maxdim = 2,
        split_order = (:up, :right, :down, :left),
    )

    @test local_density_simple(psi_reordered, center) ≈
          local_density_simple(psi_default, center) atol = 1e-10 rtol = 1e-8
    @test nearest_neighbor_density_simple(psi_reordered, center, :right) ≈
          nearest_neighbor_density_simple(psi_default, center, :right) atol = 1e-10 rtol =
        1e-8
    @test star_expectation_simple(psi_reordered, center, Hstar) ≈
          star_expectation_simple(psi_default, center, Hstar) atol = 1e-10 rtol = 1e-8
end

@testset "nontrivial D=2 fixture exercises active virtual sectors" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    Hstar = square_pxp_star_hamiltonian()
    psi_zero = _seeded_nontrivial_d2_ipeps_test(cell)

    project_star!(psi_zero, center, 0.0; projected = false, maxdim = 2)

    @test _has_active_second_virtual_sector_test(psi_zero)
    _assert_finite_star_update_state_test(psi_zero)
    @test isfinite(local_density_simple(psi_zero, center))
    @test isfinite(real(star_expectation_simple(psi_zero, center, Hstar)))

    psi_default = _seeded_nontrivial_d2_ipeps_test(cell)
    psi_reordered = _seeded_nontrivial_d2_ipeps_test(cell)
    project_star!(psi_default, center, 0.01; evolution = :imaginary, projected = true, maxdim = 2)
    project_star!(
        psi_reordered,
        center,
        0.01;
        evolution = :imaginary,
        projected = true,
        maxdim = 2,
        split_order = (:up, :right, :down, :left),
    )

    @test local_density_simple(psi_reordered, center) ≈
          local_density_simple(psi_default, center) atol = 1e-4 rtol = 1e-3
    @test nearest_neighbor_density_simple(psi_reordered, center, :right) ≈
          nearest_neighbor_density_simple(psi_default, center, :right) atol = 1e-8 rtol =
        2e-2
    @test star_expectation_simple(psi_reordered, center, Hstar) ≈
          star_expectation_simple(psi_default, center, Hstar) atol = 2e-5 rtol = 1e-1

    psi_step = _seeded_nontrivial_d2_ipeps_test(cell)
    info = project_star!(psi_step, center, 0.01; evolution = :real, projected = true, maxdim = 2)
    _assert_finite_star_update_state_test(psi_step, [info])
    @test isfinite(log_norm(psi_step))
end

@testset "D=2 update tracks normalization scale ledger" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    center = SquareCoord(5, 5)
    accumulated = 0.0

    for _ = 1:20
        info = project_star!(psi, center, 0.01; evolution = :real, projected = true, maxdim = 2)
        accumulated += sum(log, values(info.norm_factors))
        @test log_norm(psi) ≈ accumulated atol = 1e-10 rtol = 1e-10
        @test all(T -> isfinite(norm(T)), values(psi.tensors))
    end
end

@testset "repeated nontrivial D=2 star updates keep diagnostics finite" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = _seeded_nontrivial_d2_ipeps_test(cell)
    center = SquareCoord(5, 5)
    accumulated = 0.0

    @test _has_active_second_virtual_sector_test(psi)
    for _ = 1:100
        info = project_star!(psi, center, 0.005; evolution = :real, projected = true, maxdim = 2)
        accumulated += sum(log, values(info.norm_factors))
        @test log_norm(psi) ≈ accumulated atol = 1e-10 rtol = 1e-12
        _assert_finite_star_update_state_test(psi, [info])
    end

    summary = measure_simple(psi)
    _assert_finite_simple_summary_test(summary)
    @test 0 <= summary.density <= 1
    @test summary.blockade_violation >= 0
    @test summary.mean_bond_entropy >= 0
    @test summary.max_bond_entropy >= summary.mean_bond_entropy
end

@testset "repeated nontrivial D=2 split orders agree on observables" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    Hstar = square_pxp_star_hamiltonian()
    psi_default = _seeded_nontrivial_d2_ipeps_test(cell)
    psi_reordered = _seeded_nontrivial_d2_ipeps_test(cell)

    for _ = 1:25
        project_star!(
            psi_default,
            center,
            0.005;
            evolution = :real,
            projected = true,
            maxdim = 2,
            split_order = (:right, :up, :left, :down),
        )
        project_star!(
            psi_reordered,
            center,
            0.005;
            evolution = :real,
            projected = true,
            maxdim = 2,
            split_order = (:up, :right, :down, :left),
        )
    end

    _assert_finite_simple_summary_test(measure_simple(psi_default))
    _assert_finite_simple_summary_test(measure_simple(psi_reordered))
    @test local_density_simple(psi_reordered, center) ≈
          local_density_simple(psi_default, center) atol = 1e-10 rtol = 1e-8
    @test nearest_neighbor_density_simple(psi_reordered, center, :right) ≈
          nearest_neighbor_density_simple(psi_default, center, :right) atol = 1e-10 rtol =
        1e-8
    @test star_expectation_simple(psi_reordered, center, Hstar) ≈
          star_expectation_simple(psi_default, center, Hstar) atol = 1e-10 rtol = 1e-8
end

@testset "project_star accepts explicit star models" begin
    cell = PeriodicSquareUnitCell(10, 10)
    center = SquareCoord(5, 5)
    dt = 0.01

    legacy = product_square_ipeps(cell; state = :down, maxdim = 1)
    explicit = product_square_ipeps(cell; state = :down, maxdim = 1)
    project_star!(legacy, center, dt; projected = true, maxdim = 1)
    project_star!(explicit, center, dt; model = PXPStarModel(true), maxdim = 1)

    @test local_density_simple(explicit, center) ≈ local_density_simple(legacy, center) atol =
        1e-12
    @test log_norm(explicit) ≈ log_norm(legacy) atol = 1e-12

    tfim = product_square_ipeps(cell; state = :up, maxdim = 1)
    info = project_star!(tfim, center, 0.0; model = TFIMStarModel(1.0, 0.0), maxdim = 1)
    @test info.max_truncerr ≥ 0
    @test isfinite(log_norm(tfim))
end
