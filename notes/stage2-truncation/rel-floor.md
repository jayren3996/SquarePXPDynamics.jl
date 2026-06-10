# rel_floor — the conditioning knob (provisional default 1e-3)

**ADDENDUM 2026-06-10: the trajectory verdict below is LEGACY-backend-scoped.** On the
NATIVE backend (`evolve_pepskit!`) rel_floor=1e-3 is an accuracy artifact (it flattened
the D-ladder, D=3≡D=4 to 6 digits) and **rel_floor=0 is correct and stable at D=3 and
D=4**: collapse+onset tracked to ~1e-4, monotone D-ladder, no conditioning crash
(D=4 RMS 5.9e-3, best to date). See
`2026-06-10-d3-simple-update-6x6-benchmark.md` and `data/native_relfloor0_baseline.out`.
The 2026-06-02 table below (incl. the "rel_floor=0 crashed at D=4" entry) came from the
legacy `dev_relfloor_trajectory.jl` path and the since-fixed shifted-ED comparison.

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

**TRAJECTORY VERDICT (2026-06-02): 1e-3 is justified — keep it.** Re-evaluated all
floors by trajectory RMS (`scripts/dev_relfloor_trajectory.jl`), not the t=2.6
endpoint:

| floor | D=2 | D=3 | D=4 |
|---|---|---|---|
| 1e-4 | 2.87e-2 | **1.68e-2** | 1.47e-2 |
| **1e-3** | **2.62e-2** | 1.82e-2 | **9.8e-3** |
| 0 | 2.91e-2 | 2.90e-2 | crashed (conditioning) |

1e-3 wins D=2 and (clearly) D=4, is marginally behind 1e-4 at D=3, and rel_floor=0
is bad by trajectory (D=3 2.9e-2, and a high-D conditioning crash) — so a floor is
needed and adaptive-toward-0 is wrong. The 1e-4 "D=4 catastrophe" WAS
endpoint-EXAGGERATED (3.65× at t=2.6 → only 1.5× by trajectory) but 1e-3 is
genuinely better at D=4. So the default change STANDS. The floor is a conditioning
knob, not an accuracy fix; the principled fix is environment-aware truncation
(`improvement-roadmap.md`). **DONE (2026-06-02):** the revival regression gate in
`test_d_convergence.jl` now uses the trajectory RMS metric (RMS of |n(t)−n_ED(t)|
over [0,2.8], sampled every 0.2). The former t=2.6 `@test_broken err[4]<=err[3]`
endpoint inversion is replaced by PASSING `@test rms[3]<=rms[2]`, `@test rms[4]<=rms[3]`
monotonicity (measured RMS D2 2.62e-2, D3 1.82e-2, D4 9.79e-3; true gaps ~8e-3). The
endpoint was the wrong gate even though the floor choice it informed happened to be
right.
