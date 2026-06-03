using LinearAlgebra

using SquarePXPDynamics.FiniteIPEPSObservables:
    dense_state_finite,
    exact_blockade_violation_finite,
    exact_density_finite,
    exact_nearest_neighbor_expectation_finite,
    exact_one_site_expectation_finite,
    exact_pxp_energy_density_finite,
    exact_star_expectation_finite

function _finite_obs_state_after_serial_stars(nstars; step = 0.02, maxdim = 2, cutoff = 1e-12)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    for center in cell.reps[1:nstars]
        project_star!(
            psi,
            center,
            step;
            evolution = :real,
            projected = true,
            maxdim,
            cutoff,
        )
    end
    return psi
end

@testset "exact finite iPEPS observables product limits" begin
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)

    @test exact_density_finite(down) ≈ 0.0 atol = 1e-15
    @test exact_density_finite(up) ≈ 1.0 atol = 1e-15
    @test exact_blockade_violation_finite(down) ≈ 0.0 atol = 1e-15
    @test exact_blockade_violation_finite(up) ≈ 1.0 atol = 1e-15
    @test exact_pxp_energy_density_finite(down) ≈ 0.0 atol = 1e-15
end

@testset "exact finite dense state absorbs each canonical lambda once" begin
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = 2)
    for (i, key) in enumerate(sort(collect(keys(psi.link_weights)); by = k -> (k.site.y, k.site.x, String(k.dir))))
        set_link_weight!(psi, key.site, key.dir, [1.0 + i / 10, 0.25])
    end

    state = dense_state_finite(psi)
    expected = prod(link_weight(psi, key.site, key.dir)[1] for key in keys(psi.link_weights))

    @test count(x -> abs(x) > 1e-14, state) == 1
    @test state[end] ≈ expected atol = 1e-12 rtol = 1e-12
end

@testset "exact finite observables expose D2 simple measurement boundary" begin
    psi = _finite_obs_state_after_serial_stars(3)
    n = projector_up()
    nn = kron(n, n)
    zz = kron(pauli_z(), pauli_z())
    Hstar = square_pxp_star_hamiltonian()
    star_center_density = embed_one_site(n, 1, SQUARE_STAR_SITES)

    @test exact_density_finite(psi) ≈ 0.00013326224449912612 atol = 5e-15
    @test density_simple(psi) ≈ 0.000111054099003352 atol = 5e-15

    @test real(exact_one_site_expectation_finite(psi, SquareCoord(2, 1), n)) ≈
          0.0003997867121725804 atol = 5e-15
    @test local_density_simple(psi, SquareCoord(2, 1)) ≈
          0.0001999133387617944 atol = 5e-15

    @test real(exact_nearest_neighbor_expectation_finite(psi, SquareCoord(1, 1), :right, nn)) ≈
          nearest_neighbor_density_simple(psi, SquareCoord(1, 1), :right) atol = 5e-8
    @test real(exact_nearest_neighbor_expectation_finite(psi, SquareCoord(1, 1), :right, zz)) ≈
          0.9984005332366328 atol = 5e-15

    @test real(exact_star_expectation_finite(psi, SquareCoord(1, 1), star_center_density)) ≈
          0.00039994666951102603 atol = 5e-15
    @test star_expectation_simple(psi, SquareCoord(1, 1), Hstar) ≈
          exact_star_expectation_finite(psi, SquareCoord(1, 1), Hstar) atol = 5e-8
end

@testset "exact finite observables reject large cells by default" begin
    psi = product_square_ipeps(PeriodicSquareUnitCell(4, 4); state = :down, maxdim = 1)
    @test_throws ArgumentError exact_density_finite(psi)
end

@testset "boundary contractor matches dense oracle" begin
    # Product limits via the double-layer boundary path.
    cell = PeriodicSquareUnitCell(3, 3)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    up = product_square_ipeps(cell; state = :up, maxdim = 1)
    @test exact_density_finite(down; method = :boundary) ≈ 0.0 atol = 1e-12
    @test exact_density_finite(up; method = :boundary) ≈ 1.0 atol = 1e-12

    # 3x3: boundary == dense to ~1e-9 on an entangled serial-star state.
    psi3 = _finite_obs_state_after_serial_stars(3)
    @test exact_density_finite(psi3; method = :boundary) ≈
          exact_density_finite(psi3; method = :dense) atol = 1e-9

    # 4x4 (16 sites): boundary == dense across a D-ladder. The dense path is the
    # trusted reference here (D<=4 fits memory); boundary must reproduce it.
    for D in (2, 3, 4)
        cell4 = PeriodicSquareUnitCell(4, 4)
        psi = product_square_ipeps(cell4; state = :down, maxdim = 1)
        for center in cell4.reps[1:4]
            project_star!(psi, center, 0.02; evolution = :real, projected = true,
                          maxdim = D, cutoff = 1e-12)
        end
        dense = exact_density_finite(psi; max_sites = 16, method = :dense)
        bnd = exact_density_finite(psi; max_sites = 16, method = :boundary)
        @test bnd ≈ dense atol = 1e-9
        # :auto uses dense at 16 sites (the 2^N state fits); it must match.
        @test exact_density_finite(psi; max_sites = 16, method = :auto) ≈ dense atol = 1e-9
    end

    @test_throws ArgumentError exact_density_finite(psi3; method = :bogus)
end

@testset "boundary oracle handles D=6 with bounded memory (4x4)" begin
    # The boundary contractor closes physical legs locally, so its memory is
    # bounded regardless of D — it returns a finite density in [0, 1] at D=6 on a
    # 16-site torus on any host. (The dense path also reaches D=6 on a large-memory
    # host but can OOM on a constrained one; we deliberately do NOT invoke it here.)
    # Kept as a fast smoke (a few star projections, not a full revival).
    cell = PeriodicSquareUnitCell(4, 4)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    for center in cell.reps[1:6]
        project_star!(psi, center, 0.02; evolution = :real, projected = true,
                      maxdim = 6, cutoff = 1e-12)
    end
    GC.gc()
    dens = exact_density_finite(psi; max_sites = 16, method = :boundary)
    @test isfinite(dens)
    @test 0.0 <= dens <= 1.0
end
