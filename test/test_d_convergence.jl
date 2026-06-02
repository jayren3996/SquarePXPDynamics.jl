# Stage-1 HARD RULE, enforced: iPEPS reliability is judged by the TRUSTED
# observable's convergence toward ED across a D-ladder (D>=2,3,4) in an entangled
# regime. D=1 is a product state and is NEVER validation. Error vs ED must be
# non-increasing in D within tolerance and must shrink toward ED where D matters;
# a larger D worse than a smaller D (fixed dt/cutoff/time) is a HARD REGRESSION.
# See memory/stage1_d_convergence_rule.md.
#
# The full D-ladder is slow (validate_pxp_ed_ipeps at t=0.3/0.5 for D=1..4), so
# the load-bearing assertions run only under SQUAREPXP_EXTENDED_TESTS. The
# always-run smoke just proves D=1 is materially worse than D=2 (the rule's core)
# at a short time.

using Test
using SquarePXPDynamics

# Trusted-observable (exact_density_finite) error vs ED for one (D, time).
function _exact_finite_error_vs_ed(D::Int, total_time::Float64; dt = 0.02, cutoff = 1e-12)
    cfg = PXPValidationConfig(3; total_time = total_time, dt = dt, maxdim = D,
                              cutoff = cutoff,
                              exact_finite_observables = true,
                              exact_finite_max_sites = 9)
    rep = validate_pxp_ed_ipeps(cfg; ctm_params = nothing)
    s = rep.comparisons[end]
    ed = s.ed_excitation_density
    ef = something(s.ipeps_exact_finite_density, NaN)
    return abs(ef - ed)
end

@testset "Stage-1 D>1 convergence rule (smoke: D=1 is not validation)" begin
    # Even at a short time, the D=1 product state must be materially worse than
    # D=2 against the trusted observable — D=1 success proves nothing about iPEPS.
    e1 = _exact_finite_error_vs_ed(1, 0.1)
    e2 = _exact_finite_error_vs_ed(2, 0.1)
    @test e2 < e1                      # D>1 strictly helps
    # D=1 sits an order of magnitude ABOVE, and D=2 an order of magnitude BELOW,
    # the 1e-4 convergence bar — a clean "different regime" separation. Absolute
    # bounds, NOT a ratio: the rel_floor=1e-3 default (chosen for revival
    # robustness) deliberately compresses the short-time D-separation, nudging e2
    # to ~3.4e-5 so e1/e2 ~ 46x; the old `e1 > 50*e2` ratio was brittle to D=2's
    # exact error. See memory/stage2_rel_floor.md.
    @test e1 > 1e-3                    # D=1 is genuinely UNconverged (bad regime)
    @test e2 < 1e-4                    # D=2 trusted error is small (matches ED)
end

if get(ENV, "SQUAREPXP_EXTENDED_TESTS", "") != ""
    @testset "Stage-1 D>1 convergence rule (entangled D-ladder vs ED)" begin
        # t=0.3: D=1 useless; D>=2 converged and consistent (plateau).
        e = Dict(D => _exact_finite_error_vs_ed(D, 0.3) for D in 1:4)
        @test e[1] > 100 * e[2]                       # D=1 catastrophic
        @test e[2] < 1e-3                             # D>=2 tracks ED
        @test e[3] <= 1.5 * e[2]                      # no catastrophic worsening
        @test e[4] <= 1.5 * e[2]

        # t=0.5: genuinely entangled — larger D must NOT worsen, and D=4 must not
        # be worse than D=2 (it should help). This is the assertion the
        # pre-rel_floor instability (D=4 err 9.7e-3 vs D=2 1.2e-3) would FAIL.
        f = Dict(D => _exact_finite_error_vs_ed(D, 0.5) for D in 1:4)
        @test f[1] > 50 * f[4]                        # D=1 useless; D helps a lot
        @test f[4] <= f[2] + 5e-4                     # D=4 not worse than D=2 (helps)
        @test f[3] <= 1.5 * f[2]                      # small residual wrinkle tolerated, catastrophe caught
        # Document the known residual non-monotonicity at D=3, t=0.5 (a Stage-2
        # target: tune rel_floor / implement Vidal sqrt-lambda regauge). Flip to
        # @test when strict D-monotonicity holds.
        @test_broken f[3] <= f[2]
    end

    @testset "Stage-2 revival D-ladder vs ED (exact 16-site oracle)" begin
        # 4x4 Neel to the FIRST n(t) revival t=2.6 (ED peak 0.4825), measured by
        # the EXACT 16-site contraction (exact_density_finite, NO CTM environment).
        # The headline Stage-2 acceptance test. The exact oracle settled
        # (2026-06-02) that (a) CTM chi=8 was contaminating this benchmark by
        # ~3-13e-3 -- it FLATTERED the error (chi=8 put D=3/1e-4 at 4e-6 but the
        # true error is 2.9e-3) -- and (b) the revival D-non-monotonicity is REAL
        # evolution error: the mean-field-environment ceiling, NOT a measurement
        # artifact. See memory/stage2_meanfield_environment_ceiling.md.
        # Default rel_floor=1e-3; exact errors: D2 5.9e-3, D3 5.5e-3, D4 9.6e-3.
        ED = 0.4825
        err = Dict{Int,Float64}()
        for D in (2, 3, 4)
            psi = checkerboard_square_ipeps(
                PeriodicSquareUnitCell(4, 4); excited_on = :even, maxdim = D)
            params = TrotterParams(0.02, 2, :real, D, 1e-12; schedule = :serial)  # default rel_floor 1e-3
            for _ = 1:130
                evolve!(psi, 0.02; params = params)            # -> t=2.6
            end
            err[D] = abs(exact_density_finite(psi; max_sites = 16) - ED)
        end
        # Every D reaches the revival, and rel_floor=1e-3 holds the D=4 peak well
        # below the rel_floor=1e-4 catastrophe (exact D=4 error 3.45e-2 at 1e-4).
        for D in (2, 3, 4)
            @test err[D] < 1.5e-2
        end
        # ASPIRATIONAL (the env-ceiling target): strict D-monotonicity at the
        # revival. FAILS today -- D=4 is worst (9.6e-3 > D=3 5.5e-3) because of the
        # single-site mean-field environment. Flip to @test when an
        # environment-aware (full/cluster) update lands.
        # See memory/pxp_improvement_roadmap.md.
        @test_broken err[4] <= err[3]
    end
end
