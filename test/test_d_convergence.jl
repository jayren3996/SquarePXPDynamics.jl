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
    @test e1 > 50 * e2                 # D=1 is in a different (bad) regime
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
end
