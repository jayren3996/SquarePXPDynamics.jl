# D>1 validation of the PEPSKit-native backend.
#
# WHY D=1 IS INSUFFICIENT: at maxdim = 1 every star truncates back to a product
# state, so the native and legacy engines are bit-for-bit identical BY
# CONSTRUCTION and the whole truncation / SVD / leg-routing / sqrt(λ)-gauge
# machinery is bypassed. A coordinate transpose/reflection, a center↔leaf leg
# swap (invisible because all four neighbour projectors are the identical P↓), or
# a √λ-vs-λ double-count would all pass every existing D=1 native test. These
# tests exercise the native backend at D>1 against the independently ED-validated
# legacy exact oracle, plus a self-contained convention-B (source-of-truth) check.
#
# Trusted reference: `exact_density_finite(::SquareIPEPSState)` (the legacy
# CTM-free 2^N contraction, validated against ED in test_d_convergence.jl). The
# new `exact_density_finite(::PXPIPEPSState)` contracts ONLY measurement_peps.

using Test
using SquarePXPDynamics

const _n_up = ComplexF64[1 0; 0 0]   # Rydberg/excited density operator (basis 1 = up)

# Build a genuinely entangled D=2 legacy state from all-down + short serial PXP
# evolution (all-down is blockade-safe on any cell; serial works on Lx,Ly ≥ 3).
function _entangled_legacy(cell; dt = 0.1, total = 0.3, maxdim = 2)
    psi = product_square_ipeps(cell; state = :down, maxdim = maxdim)
    params = TrotterParams(dt, 2, :real, maxdim, 1e-12; schedule = :serial)
    evolve!(psi, total; params = params)
    return psi
end

@testset "native D>1 validation oracle" begin

    @testset "T0: native exact density is exact at D=1 (absolute pin)" begin
        cell = PeriodicSquareUnitCell(3, 3)
        # all-down → 0 ; x_plus → 1/2 ; checkerboard(even) → 1/2 (half excited)
        for (st, want) in ((:down, 0.0), (:x_plus, 0.5))
            nat = product_pxp_pepskit_state(cell; state = st, D = 1)
            @test exact_density_finite(nat; max_sites = 16) ≈ want atol = 1e-12
            # at D=1 the exact oracle and the mean-field simple density coincide
            @test exact_density_finite(nat; max_sites = 16) ≈ measure_simple(nat).density atol = 1e-12
        end
        nat_cb = checkerboard_pxp_pepskit_state(cell; excited_on = :even, D = 1)
        # NOTE: a 3×3 checkerboard is NOT balanced (5 excited / 9 sites), so the
        # mean density is 5/9, not 1/2 — compute it from the excitation pattern.
        want_cb = count(c -> iseven(c.x + c.y), cell.reps) / length(cell.reps)
        @test exact_density_finite(nat_cb; max_sites = 16) ≈ want_cb atol = 1e-12
    end

    @testset "T1: native exact density == legacy exact oracle at D=2 (keystone)" begin
        cell = PeriodicSquareUnitCell(3, 3)
        psi = _entangled_legacy(cell)
        # confirm it is genuinely entangled (else the test is a D=1 trap in disguise)
        @test max_bond_entropy(psi) > 1e-6
        nat = pxpipeps_from_square_ipeps(psi)
        d_leg = exact_density_finite(psi; max_sites = 16)
        d_nat = exact_density_finite(nat; max_sites = 16)
        @test d_nat ≈ d_leg atol = 1e-9
        # and the native exact density must DIFFER from the (mean-field) simple one
        # at D>1 — proving it is the gauge-invariant contraction, not the diagnostic.
        @test abs(d_nat - measure_simple(nat).density) > 1e-6
    end

    @testset "T2: convention-B — measurement_peps is the source of truth" begin
        # The exact oracle reads ONLY measurement_peps; corrupting the SUWeight
        # ledger must NOT change it, while the mean-field simple density MUST.
        cell = PeriodicSquareUnitCell(3, 3)
        nat = pxpipeps_from_square_ipeps(_entangled_legacy(cell))
        d_exact_before = exact_density_finite(nat; max_sites = 16)
        d_simple_before = measure_simple(nat).density

        corrupt = copy_state(nat)
        nontrivial = false
        for w in pepskit_weights(corrupt).data
            if length(w.data) >= 2
                w.data[2:end] .*= 2.0            # non-uniform: changes the λ² environment
                nontrivial = true
            end
        end
        @test nontrivial                          # state really has D>1 bonds to corrupt
        @test exact_density_finite(corrupt; max_sites = 16) ≈ d_exact_before atol = 1e-12
        @test abs(measure_simple(corrupt).density - d_simple_before) > 1e-6
    end

    @testset "T3: native star update vs legacy at D=2 (lossless, anisotropic)" begin
        # Anisotropic cell (Lx=3 ≠ Ly=4) so horizontal and vertical bonds are
        # genuinely distinct — a center↔leaf leg-routing swap that is invisible on
        # a C4-symmetric state changes the result here.
        cell = PeriodicSquareUnitCell(3, 4)
        psi = _entangled_legacy(cell; dt = 0.1, total = 0.2, maxdim = 2)
        nat = pxpipeps_from_square_ipeps(psi)
        # anisotropy present: some east bond weight differs from some north bond weight
        W = pepskit_weights(nat).data
        aniso = maximum(abs.(W[1, r, c].data[1] - W[2, r, c].data[1])
                        for r in axes(W, 2), c in axes(W, 3))
        @test aniso > 1e-8

        c0 = SquareCoord(2, 2)
        # one LOSSLESS star (maxdim large, cutoff 0, rel_floor 0) → no truncation,
        # so both engines apply the EXACT same gate and the exact density must match.
        psi_step = copy_state(psi)
        project_star!(psi_step, c0, 0.1; maxdim = 8, cutoff = 0.0, rel_floor = 0.0)
        nat_step = copy_state(nat)
        project_star_pepskit!(nat_step, c0; dt = 0.1, maxdim = 8, cutoff = 0.0, rel_floor = 0.0)

        @test exact_density_finite(nat_step; max_sites = 16) ≈
              exact_density_finite(psi_step; max_sites = 16) atol = 1e-9
    end

    @testset "T4: coordinate-reflection trap (per-site, vertically asymmetric)" begin
        # checkerboard on an even-Ly cell is NOT invariant under the vertical
        # reflection y → Ly−y+1, so a reflected row=Ly−y+1 ↔ row=y bug mirrors the
        # density profile and is caught site-by-site (excited 1 ↔ ground 0).
        cell = PeriodicSquareUnitCell(3, 4)
        psi = checkerboard_square_ipeps(cell; excited_on = :even, maxdim = 1)
        nat = pxpipeps_from_square_ipeps(psi)
        profile = Float64[]
        for c in cell.reps
            dl = real(exact_one_site_expectation_finite(psi, c, _n_up; max_sites = 16))
            dn = real(exact_one_site_expectation_finite(nat, c, _n_up; max_sites = 16))
            @test dn ≈ dl atol = 1e-9
            push!(profile, dn)
        end
        # the profile must be genuinely two-valued (some excited, some ground),
        # otherwise a reflection would be invisible
        @test any(>(0.9), profile)
        @test any(<(0.1), profile)
    end

    if get(ENV, "SQUAREPXP_EXTENDED_TESTS", "") != ""
        @testset "T5 (extended): native vs legacy short-time trajectory at D=2" begin
            # Full native EVOLUTION pipeline (native evolve! with scheduling +
            # multi-step + truncation) vs the legacy engine, no converter. The two
            # engines run the SAME simple-update algorithm with different SVD
            # backends, so at D=2 they keep marginally different 2-D subspaces and
            # the densities drift apart by a truncation-tie amount that ACCUMULATES
            # over the trajectory (measured: ~2.5e-7 after one step → ~5e-4 by
            # t=0.4). A gross routing / scheduling / convention bug instead shows up
            # immediately as an O(1e-2)–O(1) discrepancy at the FIRST step. So we
            # gate the first step tightly and bound the whole trajectory loosely.
            cell = PeriodicSquareUnitCell(3, 3)
            dt, maxdim = 0.05, 2
            params = TrotterParams(dt, 2, :real, maxdim, 1e-12; schedule = :serial)
            psi = product_square_ipeps(cell; state = :down, maxdim = maxdim)
            nat = product_pxp_pepskit_state(cell; state = :down, D = maxdim)
            errs = Float64[]
            for _ in 1:8
                evolve!(psi, dt; params = params)
                evolve!(nat, dt; params = params)
                d_leg = exact_density_finite(psi; max_sites = 16)
                d_nat = exact_density_finite(nat; max_sites = 16)
                push!(errs, abs(d_nat - d_leg))
            end
            @test errs[1] < 1e-5          # single-step: scheduling + gate parity (tight)
            @test maximum(errs) < 2e-3    # trajectory: tracks to truncation-tie precision
        end
    end
end
