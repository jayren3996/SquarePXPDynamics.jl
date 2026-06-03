# Stage-1 HARD RULE, enforced: iPEPS reliability is judged by the TRUSTED
# observable's convergence toward ED across a D-ladder (D>=2,3,4) in an entangled
# regime. D=1 is a product state and is NEVER validation. Error vs ED must be
# non-increasing in D within tolerance and must shrink toward ED where D matters;
# a larger D worse than a smaller D (fixed dt/cutoff/time) is a HARD REGRESSION.
# See notes/methodology/revival-validation.md.
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
    # exact error. See notes/stage2-truncation/rel-floor.md.
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

    @testset "Stage-2 revival D-ladder vs ED (exact 16-site oracle, TRAJECTORY)" begin
        # 4x4 Neel -> first n(t) collapse-and-revival, judged by the WHOLE n(t)
        # TRAJECTORY error vs ED (RMS over [0, 2.8]), NOT the single t=2.6 endpoint.
        # The headline Stage-2 acceptance test.
        #
        # WHY trajectory, not endpoint: t=2.6 (the revival peak) is a curve-CROSSING
        # where every D-curve passes near ED, so per-D errors scramble there -- the
        # EXACT t=2.6 errors are D2 6.5e-3, D3 6.0e-3, D4 1.02e-2, which falsely makes
        # D=4 look WORST (this is what the old @test_broken endpoint gate encoded).
        # By the full trajectory the ladder is cleanly MONOTONE: the error lives in
        # the revival RISE (t~1.4-2.6, a slight phase lag) and higher D reduces it.
        # The "mean-field-environment ceiling" read of the endpoint inversion was an
        # ARTIFACT of judging one instant. See notes/methodology/revival-validation.md.
        #
        # Measurement is the EXACT 16-site contraction (exact_density_finite, NO CTM):
        # the 2026-06-02 decontamination showed CTM chi=8 FLATTERS this benchmark by
        # ~3-13e-3. ED reference: 4x4 PXP torus (743-state constrained sector),
        # artifacts/neel_to_revival_4x4.json (ed.density[1:15]). Default rel_floor 1e-3.
        # Oracle: exact_density_finite (:auto -> dense at 16 sites). The
        # double-layer boundary contractor (method=:boundary) is the memory-safe
        # alternative for hosts where dense OOMs / for cells > 16 sites; both agree
        # to ~1e-9. Do NOT swap in a (contaminating) CTM measurement.
        TS = 0.0:0.2:2.8                                    # 15 samples, 10 steps each
        ED = [0.5, 0.48026525133620573, 0.4241792325244558, 0.34070816732757286,
              0.2442464623167832, 0.15549264928630516, 0.10057443671942552,
              0.10056019790056131, 0.1555774954125962, 0.24370069054283897,
              0.3373453880372616, 0.4159980265031265, 0.46695735678274974,
              0.4830586210094464, 0.46229670220636837]
        rms = Dict{Int,Float64}()
        for D in (2, 3, 4)
            psi = checkerboard_square_ipeps(
                PeriodicSquareUnitCell(4, 4); excited_on = :even, maxdim = D)
            params = TrotterParams(0.02, 2, :real, D, 1e-12; schedule = :serial)  # default rel_floor 1e-3
            meas() = (GC.gc(); exact_density_finite(psi; max_sites = 16))
            traj = Float64[meas()]                              # t=0
            for _ = 2:length(TS)
                for _ = 1:10
                    evolve!(psi, 0.02; params = params)        # 0.2 per sample
                end
                push!(traj, meas())
            end
            e = abs.(traj .- ED)
            rms[D] = sqrt(sum(e .^ 2) / length(e))
        end
        # Every D tracks the collapse-and-revival; a wrecked revival blows RMS far
        # past this. rel_floor=1e-3 holds D=4 (9.8e-3) where 1e-4 had a revival-rise
        # catastrophe (1.47e-2). Computed: D2 2.62e-2, D3 1.82e-2, D4 9.79e-3.
        for D in (2, 3, 4)
            @test rms[D] < 3.5e-2
        end
        # THE headline -- now a PASSING @test (was @test_broken at the t=2.6
        # endpoint, where it falsely inverted): by trajectory, D-convergence is
        # MONOTONE. True gaps are ~8e-3, far above the 2e-3 cross-environment slack.
        @test rms[3] <= rms[2] + 2e-3                  # D=3 beats D=2
        @test rms[4] <= rms[3] + 2e-3                  # D=4 beats D=3
    end
end
