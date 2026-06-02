# PXP dynamics-quality improvement brainstorm — synthesized roadmap (2026-06-02)

6 diverse-lens agents (environment, time-integration, ansatz, measurement,
validation, error-budget) → 32 raw ideas → synthesis. All lenses were grounded in
the four known learnings (evolution is exact-finite-correct; the ceiling is the
λ² mean-field environment; chi/regauge/floor-tuning are decisively ruled out;
PXP has an exact 2-color Trotter split), so nothing re-proposes a ruled-out fix.
Full agent output: workflow `wf_33bfb31a-411`. See
[[stage2-meanfield-environment-ceiling]], `2026-06-02-neel-to-revival-result.md`,
`2026-06-02-stage2-regauge-map.md`.

## The shape: three levers, only TWO attack the revival ceiling
- **env (the crux)** — replace the λ² mean-field environment in the truncation
  with a real one. The ONLY family that attacks the revival D-non-monotonicity.
- **measurement + validation (the scaffolding)** — make "quality" exact and
  ED-free-measurable so any env fix can be PROVEN, not just eyeballed. Highest
  tractability; the precondition for crediting the env track.
- **trotter (orthogonal)** — the exact 2-color split enables clean order-4, but
  it only cuts Trotter error, NOT the env ceiling (helps short-time, not t=2.6).

## Ranked ideas (21; rank, lever, payoff/tractability, attacks-ceiling?)
1. **[env, H/M, ✓]** Full update: reweight the bond SVD by the real CTMRG bond
   environment (build N from existing `CTMRGEnv`, minimize ‖θ−θ_trunc‖_N, regularized
   eigh + pseudo-inverse for the indefinite real-time metric) at `_split_reduced_theta:334/342`.
2. **[meas, H/H, ✗]** Make the 4×4 revival measurement EXACT (row-by-row
   boundary-MPS contraction of the 16-site double layer, bond ≤ D⁴, or raise
   `exact_finite_max_sites` 12→16). Decontaminates every test → the trusted oracle.
3. **[trotter, M/H, ✗]** Exact 2-color `:checkerboard` schedule + order-4 (Yoshida
   S4 on the verified 2-term split). Fewer truncating SVDs/dt; short-time only.
4. **[valid, H/H, ✗]** Reverse-evolution / Loschmidt echo R(t) as an ED-free
   truncation-irreversibility meter + time-to-divergence detector (scales to 5×5/6×6).
5. **[env, H/M, ✓]** Cluster / 2×2-plaquette simple update: contract the nearest
   loops EXACTLY, mean-field only on the cluster boundary. Cheaper/more stable than
   #1 (no CTM self-consistency); radius R interpolates mean-field→full update.
6. **[valid, H/H, ✗]** Per-component error budget: Richardson dt- and chi-extrap
   to split Trotter vs truncation vs measurement error at each t. The measuring stick.
7. **[env, M/M, ✓]** Belief-propagation (BP) gauge / message environment as the
   truncation metric (BP = provably loop-free-optimal local env; no CTM solve/sweep).
8. **[meas, M/H, ✗]** Env-weighted truncation-error certificate: per-step ED-free
   quality readout from the env-weighted discarded weight.
9. **[valid, H/H, ✗]** Commit the revival benchmark as a dual-regime Pareto
   regression with D-monotonicity gates + `@test_broken` peak (auto-flips green
   when a real env fix lands).
10. [meas, M/M, ✗] Return-amplitude / Loschmidt / ED-fidelity trajectory observables.
11. [env, M/L, ✓] Loop-corrected BP (leading plaquette-loop correction to the BP metric).
12. [valid, M/H, ✗] Live energy-conservation drift E(t)−E(0) tripwire (ED-free).
13. [ansatz, H/M, ✗] Constraint-resolved (Fibonacci-fusion) bond basis: encode the
    blockade structurally so forbidden near-null directions never arise.
14. [env, M/M, ✓] Imaginary-time / PSD-surrogate regularization to make the
    indefinite real-time full-update environment usable.
15. [ansatz, H/L, ✓] isoTNS / isometric-PEPS + Moses move (exact canonical form,
    true bond environment) — different ansatz CLASS, so an exact gauge EXISTS.
16. [valid, M/M, ✓] Truncation "bagging" over the exact-gate freedom (ensemble spread = error bar).
17. [env, M/L, ✓] Neural environment surrogate: learn the loop correction to λ².
18. [valid, M/M, ✗] Independent MPS/iTEBD oracle for the 4×4 cell.
19. [valid, L/M, ✗] Freeze ED golden trajectories at 5×5/6×6 (is t=2.6 peak special to 4×4?).
20. [ansatz, L/H, ✗] Larger commensurate unit cell (relieves aliasing only, not the crux).
21. [ansatz, L/L, ✗] Gaussian-fermion / parton seed ansatz.

## Ruled-out discipline held (nothing re-litigates the three negatives)
- No idea re-runs `canonicalize_simple!` expecting a dynamics fix. #7 (BP) and #15
  (isoTNS) use the word "gauge" but SURVIVE: #7 changes the DISCARDED directions
  (a pure relabel can't); #15 changes the ansatz class so an exact gauge EXISTS —
  the ingredient regauging lacked.
- No idea sweeps chi to fix the revival. #2 REMOVES CTM from 4×4; #6/#8 use chi/env
  only for error bars / readouts; #1/#14 use the env INSIDE the truncation.
- No env idea is rel_floor-tuning: #1/#5/#7/#11/#14 truncate a DIFFERENT
  (env-weighted/loop-corrected) spectrum; the "no floor is D-monotone" result
  doesn't apply to changing the metric. #13 removes directions structurally (not a floor).

## Critical gaps no lens closed (read before implementing the env track)
- **Real-time INDEFINITE metric (correctness gap, not engineering).** In real time
  the reduced environment N is non-Hermitian/indefinite (the state is not a ground
  state), so "weighted SVD" / "true-RDM truncation" is NOT a clean PSD problem.
  Only #14 confronts it (as a workaround). The right real-time truncation
  variational principle when N is indefinite is unresolved — this will bite a naive
  full-update port. Resolve the objective before coding #1.
- **No cost model.** Nobody scoped FLOPs/memory/64-core speedup for #1 vs #5 vs #7
  vs #11 vs #15. The env-track sequencing should be cost-gated; this is the single
  biggest missing input.
- **Disentangler-before-truncation** (reduce the entanglement that must be
  truncated, beyond isoTNS) — unexplored lever.
- **Is t=2.6 on 4×4 even the right target?** The whole quality def is anchored to
  one density component on one torus; finite-size/commensuration not isolated (#19/#20 peripheral).
- **Cross-observable robustness:** is D=4-worse-at-revival true for fidelity/energy
  too, or density-only? Changes the diagnosis if observable-specific.
- **Closed-loop adaptation:** the diagnostics (#4/#6/#8/#12) are passive; no
  controller adapts dt/D/env when a tripwire fires.

## Recommended next 3 (the synthesis lead's sequencing)
1. **Make the 4×4 benchmark EXACT + commit it as a regression + stand up the error
   budget (#2 + #9 + #6).** Precondition for everything: until the headline metric
   is measured with no CTM confound and encoded as a binary D-monotonicity gate,
   no env fix can be PROVEN better. Lowest-risk, highest-tractability; the
   `@test_broken` peak gate auto-flips green when a real fix lands.
   - First step: write the boundary-MPS exact contractor (bond ≤ D⁴) in
     `src/Observables.jl`, cross-check vs `exact_density_finite` on a frozen 3×3 to
     machine precision, re-run Néel-to-revival D=2,3,4 through it, compare t=2.6 to
     the chi=8 value 0.022 (match → env verdict now exact-confirmed; shift →
     quantifies prior CTM contamination).
2. **Prototype env-aware truncation reusing the existing `CTMRGEnv` in
   `_split_reduced_theta` (#1; #5 cluster as the cheaper parallel track).** The only
   lever at the crux; infra is unusually ready (`pepskit_ctmrg_context` already
   builds the env; `canonicalize_simple!` already provides the canonicalization a
   full-update step needs; hook lines documented). Let the #2/#6 oracle pick the
   winner on accuracy-per-cost. MUST resolve the indefinite-metric objective first.
   - First step: precompute the reduced bond environment N from the existing
     `CTMRGEnv` in `project_star!`, factor N = Z†Z (regularized eigh + pseudo-inverse
     on near-null dirs — mandatory, N is indefinite), env-weighted SVD of Z_L θ Z_R,
     conserve env-norm at line 342, gate behind the t=0.2 no-regression check (still
     ~1e-7 when truncerr~0).
3. **Add ED-free reverse-evolution echo R(t) + energy-drift tripwires (#4 + #12 →
   #8)** so the env prototype is validatable at 5×5/6×6 where ED is unavailable and
   the goal ultimately lives. Reuses `reverse_evolve!`/`validate_pxp_reversibility`/
   `exact_pxp_energy_density_finite`. Sensitive to exactly the irreversible-truncation
   mechanism the ceiling is made of — the property a gauge transform cannot fake.
   - First step: wire R(t)=|density(reverse_evolve!(evolve!(ψ0,t),t))−density(ψ0)| on
     4×4 Néel, small dt, D=2,3,4; confirm it spikes near t~2.6 for D=4 and →0 as
     cutoff→0 short-time; then run at 6×6 vs `artifacts/neel_ed_krylov_6.json`.
