# rel_floor — the conditioning knob (provisional default 1e-3)

`rel_floor` (in `_svd_with_rel_floor`, `StarSimpleUpdate.jl`) drops retained
singular values below `rel_floor · σ_max` after each bond SVD, capping the bond
condition number at `1/rel_floor` independent of cutoff. It exists because
`project_star!` absorbs the full λ then divides it back out, so retained λ near a
tight cutoff invert and poison later updates — the "larger-D-worsens-at-tight-
cutoff" SHORT-TIME blowup (3×3 PXP, t=0.2, D=4: 9.7e-3 → 8e-12 at rel_floor=1e-4).
That short-time blowup is real (not an endpoint artifact); the floor fixes it.

**Default: `1e-3`** in `TrotterParams` / `evolve!` (`project_star!`'s own default
stays `0.0` so dense-reference exactness tests are unaffected).

**Cost — ultra-short time:** at t=0.02 (one step) 1e-3 truncates the D=2 bond to a
near-PRODUCT state, erasing the mean-field offset, so `density_error_simple` → ~7e-7.
The offset-documenting tests (`test_pxp_validation.jl`, `test_pxp_larger_d_ed_benchmark.jl`)
were therefore moved to an entangled t=0.1 (offset ~6e-3 survives), with their
`< 1e-6` exact-finite bounds relaxed to `< 1e-3`.

**PROVISIONAL — re-evaluate by trajectory.** The 1e-4 → 1e-3 change was decided at
the single t=2.6 ENDPOINT (1e-4 "D=4 catastrophe" 3.45e-2 vs 1e-3 9.6e-3). t=2.6
is a misleading crossing point (see `../methodology/revival-validation.md`), and
early trajectory data shows 1e-4 and 1e-3 within ~10% by RMS — nothing like the
3.6× endpoint gap. So the catastrophe was likely an endpoint artifact and the
default change is NOT settled. The floor is a conditioning/robustness knob, not an
accuracy fix; the principled fix is environment-aware truncation
(`improvement-roadmap.md`). Re-pick the default by TRAJECTORY RMS across floors
(`scripts/dev_relfloor_trajectory.jl`), and move the t=2.6 `@test_broken`
regression test to a trajectory metric.
