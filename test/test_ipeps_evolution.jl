using ITensors
using LinearAlgebra

function _seeded_nontrivial_d2_evolution_ipeps_test(cell)
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

function _assert_finite_simple_summary_evolution_test(summary)
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

@testset "Trotter parameter validation" begin
    @test_throws ArgumentError TrotterParams(0.0, 1, :real, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(Inf, 1, :real, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(NaN, 1, :real, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 3, :real, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :bad, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, 0, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, 1, -1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, 1, Inf)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, 1, NaN)
    @test_throws ArgumentError TrotterParams(
        0.1,
        1,
        :real,
        1,
        1e-12;
        schedule = :bad,
    )
    @test_throws ArgumentError TrotterParams(
        0.1,
        1,
        :real,
        1,
        1e-12,
        (:right, :up, :left, :left),
    )
    @test_throws ArgumentError TrotterParams(0.0, 1, :real, true, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 3, :real, true, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :bad, true, 1, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, true, 0, 1e-12)
    @test_throws ArgumentError TrotterParams(0.1, 1, :real, true, 1, -1e-12)
end

@testset "Trotter schedules" begin
    p1 = TrotterParams(0.1, 1, :real, 1, 1e-12)
    @test trotter_sequence(p1) == [(1, 0.1), (2, 0.1), (3, 0.1), (4, 0.1), (5, 0.1)]

    p2 = TrotterParams(0.1, 2, :real, 1, 1e-12)
    @test trotter_sequence(p2) == [
        (1, 0.05),
        (2, 0.05),
        (3, 0.05),
        (4, 0.05),
        (5, 0.1),
        (4, 0.05),
        (3, 0.05),
        (2, 0.05),
        (1, 0.05),
    ]
end

@testset "zero-time evolution does not mutate state" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
    weights_before = deepcopy(psi.link_weights)
    density_before = density_simple(psi)
    blockade_before = blockade_violation_simple(psi)
    params = TrotterParams(0.1, 2, :real, 1, 1e-12)

    log = evolve!(psi, 0.0; params = params)

    @test log.nsteps == 0
    @test isempty(log.layer_infos)
    @test log.log_norm_before == 0.0
    @test log.log_norm_after == 0.0
    @test log.log_norm_delta == 0.0
    @test density_simple(psi) ≈ density_before
    @test blockade_violation_simple(psi) ≈ blockade_before
    @test psi.link_weights == weights_before
end

@testset "evolution rejects non-integer step counts" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.1, 1, :real, 1, 1e-12)

    @test_throws ArgumentError evolve!(psi, 0.25; params = params)
    @test_throws ArgumentError evolve!(psi, -0.1; params = params)
    @test_throws ArgumentError evolve!(psi, Inf; params = params)
    @test_throws ArgumentError evolve!(psi, NaN; params = params)
end

@testset "evolution requires five-color-compatible unit cells" begin
    psi_bad = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    params = TrotterParams(0.1, 1, :real, 1, 1e-12)

    @test_throws ArgumentError evolve!(psi_bad, 0.1; params = params)
end

@testset "serial sweep evolution accepts non-five-color-compatible unit cells" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.02, 1, :real, 1, 1e-12; schedule = :serial)

    log = evolve!(psi, 0.02; params = params)

    @test log.params.schedule === :serial
    @test log.nsteps == 1
    @test length(log.layer_infos) == length(cell.reps)
    @test all(length(layer) == 1 for layer in log.layer_infos)
    @test [only(layer).center for layer in log.layer_infos] == cell.reps
    @test isfinite(log.max_truncerr)
    @test all(isfinite, values(all_bond_entropies(psi)))
end

@testset "second-order serial sweep reverses center order" begin
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.02, 2, :real, 1, 1e-12; schedule = :serial)

    log = evolve!(psi, 0.02; params = params)

    expected_centers = [cell.reps; reverse(cell.reps[1:(end - 1)])]
    @test length(log.layer_infos) == length(expected_centers)
    @test [only(layer).center for layer in log.layer_infos] == expected_centers
    @test isfinite(log.max_truncerr)
end

@testset "one first-order step from all-down state" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.02, 1, :real, 1, 1e-12)

    log = evolve!(psi, 0.02; params = params)

    @test log.nsteps == 1
    @test length(log.layer_infos) == 5
    @test log.log_norm_before == 0.0
    @test log.log_norm_after ≈ log_norm(psi) atol = 1e-12
    @test log.log_norm_delta ≈ log.log_norm_after - log.log_norm_before atol = 1e-12
    @test all(
        length(layer) == length(update_centers(cell, color)) for
        (layer, color) in zip(log.layer_infos, 1:5)
    )
    @test isfinite(log.max_truncerr)
    @test log.max_truncerr >= 0
    @test isfinite(log.max_bond_entropy)
    @test isfinite(log.mean_bond_entropy)
end

@testset "one second-order step from all-down state" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.02, 2, :real, 1, 1e-12)

    log = evolve!(psi, 0.02; params = params)

    @test log.nsteps == 1
    @test length(log.layer_infos) == 9
    @test isfinite(log.max_truncerr)
    @test all(isfinite, [log.max_bond_entropy, log.mean_bond_entropy])
end

@testset "repeated small steps remain finite" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.01, 1, :real, 1, 1e-12)

    log = evolve!(psi, 0.05; params = params)

    @test log.nsteps == 5
    @test isfinite(log.max_truncerr)
    @test all(isfinite, values(all_bond_entropies(psi)))
    for lambda in values(psi.link_weights)
        @test norm(lambda) ≈ 1 atol = 1e-10
    end
end

@testset "imaginary-time smoke test" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = TrotterParams(0.01, 1, :imaginary, 1, 1e-12)

    log = evolve!(psi, 0.01; params = params)

    @test log.nsteps == 1
    @test isfinite(log.max_truncerr)
    @test all(isfinite, values(all_bond_entropies(psi)))
end

@testset "D=2 evolution smoke test" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    params = TrotterParams(0.01, 1, :real, 2, 1e-12)

    log = evolve!(psi, 0.01; params = params)

    @test log.nsteps == 1
    @test isfinite(log.max_truncerr)
    @test all(length(lambda) <= 2 for lambda in values(psi.link_weights))
    @test all(
        isapprox(norm(lambda), 1; atol = 1e-10) for lambda in values(psi.link_weights)
    )
end

@testset "nontrivial D=2 evolution tracks split normalization ledger" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = _seeded_nontrivial_d2_evolution_ipeps_test(cell)
    params = TrotterParams(0.005, 1, :real, 2, 1e-12)

    evolution_log = evolve!(psi, 0.005; params = params)

    accumulated = sum(
        sum(Base.log, values(info.norm_factors)) for
        info in Iterators.flatten(evolution_log.layer_infos)
    )
    @test evolution_log.nsteps == 1
    @test length(evolution_log.layer_infos) == 5
    @test evolution_log.log_norm_delta ≈ accumulated atol = 1e-10 rtol = 1e-12
    @test evolution_log.log_norm_after ≈ log_norm(psi) atol = 1e-12
    @test all(T -> isfinite(norm(T)), values(psi.tensors))
    @test all(lambda -> all(isfinite, lambda), values(psi.link_weights))
    @test all(
        lambda -> isapprox(norm(lambda), 1; atol = 1e-10),
        values(psi.link_weights),
    )
    _assert_finite_simple_summary_evolution_test(measure_simple(psi))
end

@testset "evolution convenience constructor delegates to parameterized method" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)

    log = evolve!(psi, 0.01; dt = 0.01, order = 1, evolution = :real, projected = true)

    @test log.params == TrotterParams(0.01, 1, :real, 1, 1e-12)
    @test log.nsteps == 1
end

@testset "evolve accepts explicit static model protocol" begin
    cell = PeriodicSquareUnitCell(10, 10)
    params = TrotterParams(0.01, 1, :real, 1, 1e-12)

    legacy = product_square_ipeps(cell; state = :down, maxdim = 1)
    explicit = product_square_ipeps(cell; state = :down, maxdim = 1)

    legacy_log = evolve!(legacy, 0.01; dt = 0.01, order = 1, evolution = :real, projected = true)
    explicit_log = evolve!(
        explicit,
        0.01;
        params = params,
        protocol = StaticModel(PXPStarModel(true)),
    )

    @test explicit_log.nsteps == legacy_log.nsteps
    @test explicit_log.max_truncerr ≈ legacy_log.max_truncerr atol = 1e-12
    @test log_norm(explicit) ≈ log_norm(legacy) atol = 1e-12
end

@testset "legacy TrotterParams constructor remains accepted" begin
    old = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    current = TrotterParams(0.01, 1, :real, 1, 1e-12)
    @test trotter_sequence(old) == trotter_sequence(current)
    @test old.dt == current.dt
    @test old.order == current.order
    @test old.evolution == current.evolution
    @test old.projected === true

    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    log = evolve!(psi, 0.01; params = old)
    @test log.params == current
    @test log.nsteps == 1
    @test isfinite(log.max_truncerr)

    explicit = product_square_ipeps(cell; state = :down, maxdim = 1)
    @test_throws ArgumentError evolve!(
        explicit,
        0.01;
        params = old,
        protocol = StaticModel(PXPStarModel(true)),
    )

    old_unprojected = TrotterParams(0.01, 1, :real, false, 1, 1e-12)
    @test old_unprojected.projected === false
    legacy_unprojected = product_square_ipeps(cell; state = :down, maxdim = 1)
    explicit_unprojected = product_square_ipeps(cell; state = :down, maxdim = 1)
    legacy_unprojected_log = evolve!(legacy_unprojected, 0.01; params = old_unprojected)
    explicit_unprojected_log = evolve!(
        explicit_unprojected,
        0.01;
        params = current,
        protocol = StaticModel(PXPStarModel(false)),
    )
    @test legacy_unprojected_log.max_truncerr ≈ explicit_unprojected_log.max_truncerr atol =
        1e-12
    @test log_norm(legacy_unprojected) ≈ log_norm(explicit_unprojected) atol = 1e-12
end
