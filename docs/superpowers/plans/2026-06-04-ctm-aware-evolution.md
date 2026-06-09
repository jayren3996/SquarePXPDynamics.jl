# CTM-Aware (Full-Environment) Evolution of the Five-Site Square-PXP iPEPS — Stage C/D Deepening

> Status: design plan (2026-06-04). Deepens Stage C/D of
> `docs/superpowers/specs/2026-06-04-pepskit-native-ctm-aware-migration-design.md`.
> Produced by a 19-agent physics-first workflow (2 recon → 6 angle analyses →
> 9 adversarial verifiers + 1 completeness critic → synthesis), then independently
> spot-checked against PEPSKit 0.7.0 source (see the Addendum at the end).
> This is a PLAN, not implementation — no production code has changed. The first
> actionable step is **PR-0 (gate G0)**, a zero-infrastructure discriminating
> experiment that decides the entire staging before any environment code is written.

## 1. Physical thesis

The native star simple update is ED-near-exact through the collapse (t<1.4) at every D, then accumulates a *structural* phase lag on the revival rise (t≈1.4–2.6) that shrinks monotonically with D (4×4 Néel trajectory RMS over [0,2.8]: D2 2.62e-2 → D3 1.82e-2 → D4 9.79e-3, recorded in `notes/methodology/revival-validation.md:18-22` and enforced by `test/test_d_convergence.jl`). The lag is not a truncation-metric error: it is the price of applying every Trotter gate against the single-site λ² Bethe/tree environment (`PEPSKitStarUpdate.jl:280` sqrt(λ) leaf absorb, `:290` center contracted bare, `:299` gate applied in that environment). **The single sharpest physics reason CTM-aware evolution is the right fix:** the revival is the rebuilding of Néel order, which is a *closed-loop* (plaquette) correlation on the square lattice, and a tree environment severs every loop by construction — so the deficit is missing representable correlation in the *environment the gate sees*, not in the bond it truncates. The fix must restore loops to the environment **during gate application**, not merely reweight the post-gate SVD.

## 2. Physics diagnosis — what λ² discards, and the Stage-C-vs-D verdict

**What the λ² Bethe environment discards.** Each leaf is dressed by `absorb_weight(..., inv=false)` = sqrt(λ) on its three outward bonds (`PEPSKitStarUpdate.jl:280`), the center is contracted with bare reduced cores (`:294-297`), and the gate is applied to the full 5-physical-leg object (`:299`). The environment seen by the cross is therefore a *product* of independent arm spectra λ²: every path that would close a loop around the central plaquette is replaced by an uncorrelated mean field. That is exactly the antiferromagnetic loop response that rebuilds Néel order on the rise.

**The Stage-C-vs-D verdict (adversarial pass: all three lenses REFUTE bare Stage C).** The truncation hook `_split_bond` (`PEPSKitStarUpdate.jl:207-232`) is a bare `tsvd(t)` (identity/Frobenius metric, no environment `N`), and it reports the relative discarded weight ε at `:225`. A *reweighted* SVD can only redistribute that ε. The anatomy probe (`/tmp/probe2.jl`, `/tmp/probe4.jl` — **ephemeral, see §5/§6 gate G0**) measured the density error at 95–350× the discarded weight throughout the rise at both D=2 and D=4 (e.g. t=2.0: D2 |err|=5.6e-2 vs ε≈1.6e-4; D4 |err|=1.14e-2 vs ε≈1.2e-4). With ~1e-4 of weight to move against a ~1e-2 error, the kept D-subspace is already near-optimal under any metric. Corroborated by the chi 8→16 control: at fixed D=4 the revival error moved 0.0221→0.0241 because only *measurement* changed, not *evolution* (`notes/stage2-truncation/2026-06-02-neel-to-revival-result.md:38-44`). And regauging (`canonicalize_simple!`) is a pure gauge transform, empirically a no-op for dynamics (`notes/stage2-truncation/2026-06-02-stage2-regauge-map.md:71-101`).

**Verdict:** *A bare bond-environment-reweighted star-SVD (Stage C as a metric swap on `_split_bond`) CANNOT close the revival lag, and shipping it would be misread as "CTM-aware truncation does not help."* The error enters via the environment the gate is applied in, and is multi-site/loop-correlation by nature. Two consequences for staging, both departures from the original spec's "Stage C = reweighted SVD" reading:

1. The minimum *viable* environment fix is **"Stage C done right"**: re-apply the gate against the real (loop-carrying) environment and optimize the truncated bond against that environment — i.e. change the environment the gate sees, not the discard ranking. Even this is a *per-arm* (sequential, uncoupled-cross) approximation and is predicted to **plateau above the simple-update D-ladder**.
2. The principled fix that captures the joint plaquette loop is **Stage D** (5-tensor cross patch optimized jointly) or a **2×2/3×3 cluster** update. These are the only families that keep the four cross arms coupled.

**The cheap discriminating experiment to run FIRST (no new infrastructure).** Before any environment code is written, settle C-vs-D quantitatively on the 4×4 Néel cell using the exact oracle `exact_density_finite` / `_native_dense_state` (`PEPSKitObservables.jl:474-553`, confirmed present, ≤16 sites, `:dense`):

> **G0 — Frozen-state exact-best-rank-D diagnostic.** At mid-rise (freeze the native D-state at a few t∈[1.6,2.2]), build the exact 2^16 statevector, apply ONE exact square-star gate to it, and compute the exact best rank-D truncation of the gated state across each center–leaf bond by brute-force SVD of the exact reduced bond. Compare three densities: (i) exact untruncated post-gate, (ii) exact-best rank-D, (iii) the simple-update result.
> - If (ii)≈(i) but (iii) lags both → the error is the EVOLUTION/ENVIRONMENT; bare Stage C is dead; go to Stage-C-done-right / Stage D / cluster.
> - If (ii) materially beats (iii) → a better truncation metric has real headroom; bare Stage C is worth a milestone.

The §1 evidence predicts the first branch. G0 is ~a dozen exact contractions and reuses the already-validated oracle, so it is the first PR (§6).

## 3. The real-time indefinite-metric problem and the CHOSEN scheme

**The contradiction, resolved.** The migration spec calls the real-time metric `N` "non-Hermitian/indefinite, ill-posed — a correctness gap, not engineering" (design doc :177-187). The roadmap states the opposite: "the norm environment ⟨ψ|ψ⟩ is PSD even in real time" (`notes/stage2-truncation/improvement-roadmap.md:66-67`). **The code settles it and the spec commits a category error.** What PEPSKit's truncators actually minimize is the *norm/fidelity* object: `bondenv_fu` builds `N` from the bra-ket double layer `X'X`, `Y'Y` (`benv_ctm.jl`), and `cost_function_als` minimizes ‖|ab⟩ − |a'b'⟩‖²_N (`als_solve.jl`). The reduced environment of ⟨ψ|ψ⟩ is a sum of |amplitude|² and is **PSD for any wavefunction, real time or imaginary**. Minimizing ‖θ − θ_trunc‖²_N is a genuine bounded-below least-squares problem; the full complex amplitudes (and thus the revival interference phases) are retained — the "clamping kills interference" worry applies to the *energy/effective-H* environment, which **no scheme here uses**. *This must be written into the spec/roadmap as a correction before coding (§5 P0).*

**The residual indefiniteness is finite-χ noise only.** A CTM-derived `N` is PSD only as χ→∞; at finite χ it picks up small negative eigenvalues, which `positive_approx` (`gaugefix.jl:14`, confirmed: eigh of (N+N')/2, global sign flip, clamp-negatives, sqrt) regularizes. This is a small numerical clamp, not a physics compromise — **provided** the clamped negative weight is monitored (§5 gate, §7 self-deception guard).

**CHOSEN scheme: NTU-style PSD finite-cluster environment as the primary engine; CTM full-environment as the large-χ accuracy ceiling/diagnostic only.**
- **Primary (well-posed by construction):** build the bond/star environment as the *exact* bra-ket overlap `N = ⟨patch|patch⟩` of a finite closed neighborhood around the cross (the 5 star sites plus a one-plaquette ring), contracted with **no** boundary/CTM fixed point. `N` is then Hermitian-PSD *term-by-term* in real time exactly as in imaginary time — `⟨δ|N|δ⟩ = ‖patch_with_δ‖² ≥ 0` by construction, with **no positivize-by-fiat step**. This is the only scheme whose well-posedness is structural rather than asserted, and it sidesteps the entire CTM-trust apparatus (residual mislabel, χ-flattering, staleness).
- **Why CTM-FU is demoted (not rejected):** it is the only path to correlations beyond the patch radius, but it (a) needs `positive_approx` clamping (acceptable but must be monitored), (b) requires a verified-converged large-χ environment the project deliberately routes around (χ=8 flatters the revival by ~the signal size), and (c) needs a real convergence judge that does not exist yet (§5 P1). Keep it as a chi-ladder accuracy ceiling and as the >patch-radius correlation probe — never as the per-step truncation metric until P1 lands.
- **Parallel fallback (cheapest well-posed):** the 2×2/3×3 cluster update inherits PSD-ness from the existing tree boundary (exact interior, mean-field edge), needs no environment metric inversion at all, and is the lowest-risk first multi-site test of the loop-correlation hypothesis. Cost-gate NTU against it; if it matches NTU on accuracy-per-cost, it is the winner.

**Rejected outright:** (a) bare reweighted bond SVD (§2 verdict — measurably fails); (b) truncating against the *indefinite energy* environment (cost function unbounded below — no regularization makes it a minimization); (c) regularized-eigh/pseudo-inverse on a symmetrized indefinite `N` as the *primary* engine (it is the symptom, not a cure — only needed for the demoted CTM-FU path's finite-χ noise).

## 4. Environment construction for the five-site star, in PEPSKit 0.7.0 terms

### 4a. Stage-C-done-right (per-arm full-update bond step) — where it hooks

The gate is applied to the full 5-site θ *first* (`PEPSKitStarUpdate.jl:299`), so the per-arm peel does **not** collapse the cross at the gate; only the post-gate truncation is per-bond. Replace the bare `_split_bond` call at `PEPSKitStarUpdate.jl:322` with an environment-weighted bond truncation:

1. Reuse the per-arm QR-isometry `Q` already produced at `:284-285` (`leftorth(Tp)`) as the leaf-side `X`/`Y` that `bondenv_fu`/`_qr_bond` expect (PEPSKit `PEPSOrth` rank-4).
2. **NTU variant (primary):** build `N_dir` as the exact bra-ket contraction of the cross + one-ring patch with the active center–leaf bond left open (reuse PEPSKit's `X'X`/`Y'Y` double-layer stacking from `bondenv_fu`, but substitute a *finite patch* for the CTM corners/edges). PSD by construction — no `positive_approx`.
   **CTM variant (ceiling only):** `benv = bondenv_fu(r, c, X, Y, env)` → `Z = positive_approx(benv)` → optional `fixgauge_benv(Z, a, b)` to crush cond(N).
3. Truncate with the shipped `bond_truncate(a, b, N, ALSTruncation(; trunc = truncrank(D) & truncerror(; atol)))` (`bond_truncation.jl`; `ALSTruncation`/`FullEnvTruncation` are the only **exported** truncators).
4. Write back the returned `s` as the SUWeight slot via the existing symmetric `absorb_s` sqrt-split idiom (`PEPSKitStarUpdate.jl:230`), preserving convention B (§5 P3).

*Note the vertical arms:* `bondenv_fu` is a horizontal half-infinite band; up/down arms need a `rotl90`'d environment or a column-band analog (open item, §8). *Note staleness:* the four arms peel sequentially and each peel changes the center, so arm k+1's environment is stale vs arm k — rebuild `N_dir` per arm from the same patch/env (acceptable to leading order), and quantify against G0.

**Replace deprecated `_bond_trunc`** (`PEPSKitStarUpdate.jl:196-202` uses `truncdim`/`truncbelow`, deprecated aliases) with `truncrank`/`trunctol` so it composes with the `bond_truncate` strategies via `&`.

**Cost (per bond):** env build dominated by `chi²·D⁶`–`chi²·D⁸` (CTM) or the finite-patch double layer (NTU, `~D^{2z}`, z = ring connectivity); ALS inner solve `~D⁶` (`KrylovKit.linsolve` on a D²×D² normal-equation operator). At D=2–4, chi≤32 this is in-core (≤~10⁷ working set); the CTM-FU cost is dominated by the `leading_boundary` solve multiplicity, not the bond solve.

### 4b. Stage D — full 5-tensor cross patch ALS

Optimize all five updated cores (new center C′ + 4 reduced leaves) **jointly** against the closed plus-star environment `E₊`, applying the gate via the exact M7 MPO so the patch contraction never materializes the dense rank-10 gate.

- **M7 gate factorization (verified exact, rank-2-per-cut).** `star_gate_mpo` (`PEPSKitStarSnake.jl:78-106`) builds the cross gate by successive SVD along `(center,right,up,left,down)`, keeping only numerically-nonzero singular values (`tol = max(s)·size·eps·10`, `:90`), so `reconstruct_star_gate` is exact and `bond_dims` reports the true per-cut Schmidt rank. **EXPERIMENTAL, not wired in** (`:6-15`). The rank-2 claim and snake **path-independence** both lack a committed contract test — both are hard prerequisites before the MPO is used to *apply* (not just store) the gate (§8 open item, §7 G6).
- **`E₊` has NO PEPSKit primitive.** Author a bespoke `@autoopt @tensor` over the 8 cells surrounding the plus, modeled on `bondenv_fu`'s `@tensor` (`benv_ctm.jl`), leaving the 4 outward leaf bonds (bra+ket) open. For the NTU variant this is a finite closed patch (PSD); for the CTM variant it stitches 4 `EnlargedCorner`s (`sparse_environments.jl`) + 4 edge tensors. **Aliasing hazard:** the plus is a 3×3 footprint; on the 4×4 validation cell a CTM `E₊` with a neighbor ring spans ≥5 cells per axis and may self-alias under the mod1 1-cell-renormalized `CTMRGEnv` — *unverified geometric blocker* (§5 P5, §8). The NTU finite patch does not have this problem (it is not periodically renormalized).
- **ALS sweep:** hold 4 cores fixed, solve the 5th from its normal equation `R x = S` (`R` = `E₊` closed on the other 4 ket+bra cores, `S` = `E₊` closed on the M7-gated target); cycle center→R→U→L→D until fidelity converges. PEPSKit's `_solve_ab`/`linsolve`/`cost_function_als` (`als_solve.jl`) are reusable per single-core subproblem by reshaping to the rank `_solve_ab` expects (the shipped machinery is hard-wired to {2,1}/{1,2} bond tensors — per-core reshape needed, possibly a bespoke 5-core solver; open item §8).
- **Must not degenerate to 4 uncoupled bond solves** — that silently loses the multi-site correlation Stage D exists to capture (contract-test that the joint optimum differs measurably from per-bond).
- **Cost:** `E₊` close `~chi²·D⁶` (CTM) with the M7 rank-2 slicing keeping the gated-target intermediate at `~2·D²` instead of `32·D²`; per-core solve `~D⁸` (center), `~D⁶` (leaf). In-core at D≤4, chi≤32 (~1 MB transient `E₊`), but ~10–50× the Stage-C cost and **uncosted at D=4** — §5 P4 cost model is a hard gate.

## 5. Hard PREREQUISITES (each with a pass/fail gate)

| ID | Prerequisite | Pass/fail gate |
|---|---|---|
| **P0** | **Resolve PSD-vs-indefinite in the docs.** Norm env ⟨ψ\|ψ⟩ is PSD in real time; the spec's "indefinite N ill-posed" framing conflates it with the unused energy env. | Spec design doc :177-187 corrected (or roadmap :66-67 corrected) and a one-line note added: the truncation objective is the PSD norm-fidelity, `cost_function_als`/`bondenv_fu`. **No environment code merges until one doc is fixed.** |
| **P1** | **Real CTM residual η (CTM-FU path only).** `leading_boundary` info has NO η/converged/iteration field — only `truncation_error`, `condition_number`. The legacy adapter aliases `truncation_error` into the `residual` slot (`PEPSKitMeasurements.jl:493`, confirmed) and gates `accepted` on it (`:497-499`). The native context stores info raw and never inspects convergence (`PEPSKitObservables.jl:181,199-210`). | Native diagnostics compute `η = calc_convergence(env_new, env_old)` (`ctmrg.jl:191-207`, `@non_differentiable`, confirmed) across iterations; store `η`, `truncation_error`, `condition_number` as **three distinct fields**; gate `converged := η ≤ tol`; never let `truncation_error` reach `residual`. Contract test: post-hoc η matches the loop's final η on a converged env AND reports large η on a deliberately under-converged (`maxiter=2`) env. |
| **P2** | **χ-sufficiency margin (CTM-FU path only).** χ=8 flatters the revival by ~the signal size. | χ-ladder ({8,16,24,32} or {kD²}) at fixed D; require the **CTM-free** trajectory RMS (`exact_density_finite`) to be monotone-decreasing AND plateau, with η≤tol at the plateau, and env chi-drift ≤ 0.1·ε_trunc(D) (ε from `PEPSKitStarUpdate.jl:225`). Require monotone-decreasing drift across ≥3 χ points before trusting the margin. |
| **P3** | **Convention-B / ledger reconciliation.** Post-environment-weighted truncation, the bond weight is NOT a Schmidt spectrum. `measure_simple`'s `_site_array` double-uses the ledger on top of convention-B Γ; `measure_simple` is mean-field only. | After `bond_truncate` returns optimized a,b, do a bare local SVD a·b=USV purely to (i) write back symmetrically (`U·sqrt(S)`, `sqrt(S)·V`, mirroring `absorb_s`) so `state.peps` stays convention-B, and (ii) store S in the SUWeight slot as a **diagnostic-only** local bond spectrum (relabel project-wide). Demote `measure_simple` to "tree diagnostic, exact only for bare-SVD states" and **exclude it from every gate**. Forbid (or explicitly test) mixed bare-SU/CTM-aware schedules (the deabsorb floor `eps^(3/4)·λmax` divides by a non-Schmidt λ — legacy near-zero-division hazard reborn). Gate: convention-B equivalence test re-run after env-aware truncation passes. |
| **P4** | **FLOPs/memory cost model at D=2,3,4 × χ** (the brainstorm's "single biggest missing input"). | A written model pricing the CTM-solve multiplicity (per-colour × χ-ladder × D-ladder × nsteps) AND the dt-vs-staleness coupling (does holding the env across the 9 order-2 substeps degrade order-2 Trotter cancellation to order-1? — does it force per-colour rebuild = 5× solves or smaller dt?). Stage D blocked until the model shows the patch ALS is bounded at D=4 at the revival peak. |
| **P5** | **`E₊` geometry on the validation cell** (Stage D only). | Construct a plus-star env on a 4×4 `CTMRGEnv` and check for index self-aliasing (`bondenv_fu` already spans 4 columns `cm1..cp2`). If it aliases, the NTU finite patch must be used on 4×4 (CTM `E₊` validated only on ≥6×6 where no native exact oracle exists). Contract-test `E₊` against brute-force `ncon` at D=1,2. |
| **P6** | **Blockade-constraint defense under env truncation** (new, from completeness critic). The gate is `P·U` with LEFT projection (`PXPModel.jl:200-211`); `:sandwich` (P·U·P) exists. `[P,U]=0` holds on a *clean* constrained state — an env-weighted truncation re-ranks by a non-product metric and can leak blockade-violating amplitude the left-only P does not re-clean. No D>1 trajectory blockade gate exists. | Decide P·U vs P·U·P for env-aware steps; add a D>1 **trajectory** blockade-violation tripwire (cheap, gauge-invariant via `blockade_violation_ctm_pepskit`, currently unused as a tripwire). Gate: blockade violation stays ≤ the bare-SU baseline across [0,2.8]. |

## 6. Staged incremental plan (TDD, smallest reviewable PRs, legacy retained)

Each stage: **red gate → minimal impl → audit note in `notes/stage2-truncation/`**. The legacy ITensors backend and the native bare-SVD `evolve_pepskit!` are retained unchanged throughout; env-aware truncation is opt-in behind a new control bundle on `TrotterParamsPK` (the driver `PEPSKitEvolution.jl:143-189` currently has no env hook).

- **PR-0 (FIRST, small, no production code change): the G0 discriminating experiment.** A `scripts/`/test harness that runs §2's frozen-state exact-best-rank-D diagnostic at D=2,4 on the 4×4 cell at t∈{1.6,1.8,2.0,2.2}, using only `exact_density_finite`/`_native_dense_state` + the exact star gate. **Red gate:** assert the three densities and report (ii)−(i) and (iii)−(i). This is the decision gate for the entire staging: it either kills bare Stage C or revives it. Also persist the ephemeral `/tmp/probe*.jl` error/discard-ratio numbers into a committed artifact (they are the load-bearing evidence and currently exist only in `/tmp`).
- **PR-1 (CTM-trust prerequisite, P1): real η in the native context.** Red gate = the P1 contract test (η matches loop on converged env, large on `maxiter=2`). Minimal impl = `calc_convergence` shim + native `CTMRGDiagnostics` with three distinct fields. **Rewrite the tests that assert the mislabel** (`test/test_pepskit_measurements.jl:167,340` assert `diag.residual == truncation_error`; `test/test_ctm_trust.jl:286-318` `:missing_residual`/`:residual_too_large`) to assert η semantics — otherwise the bug is frozen by the suite.
- **PR-2 (P3 + shims): convention-B write-back + private-primitive shims.** Shim + contract test each PEPSKit private used (`positive_approx`, `bondenv_fu`, `bond_truncate`, `_qr_bond`, `absorb_s`, `flip_svd`, `calc_convergence`); pin PEPSKit 0.7.0. Local-SVD diagnostic-ledger re-extraction; relabel `measure_simple`; convention-B equivalence test under a non-Schmidt write-back. Add the P6 trajectory blockade tripwire.
- **PR-3 (first env-aware truncation, decided by PR-0): NTU per-arm bond step.** Only if PR-0 leaves headroom for a bond-level fix OR as the well-posed minimum-viable env step. Build the finite-patch bra-ket `N_dir`, hook at `PEPSKitStarUpdate.jl:322`, write back via PR-2. **Red gate:** the §7 G1 fixed-D lag-reduction + G2 chi(/radius)-monotonicity, measured by the CTM-free oracle. Expect a plateau above the simple-update D-ladder (per the §2 verdict) — that plateau is itself the quantitative motivation for Stage D.
- **PR-4 (parallel track): 2×2/3×3 cluster update.** Bespoke cluster contraction (no env inversion, PSD by inheritance). Same G1/G2 gates. Cost-gate against PR-3.
- **PR-5 (Stage D, gated on P4+P5 and a PR-3/PR-4 plateau): M7 wiring + 5-tensor patch ALS.** First a sub-PR wiring/contract-testing the M7 MPO application (rank-2 + path-independence, §7 G6). Then `E₊` (contract-tested vs `ncon`, P5) + the joint ALS. **Red gate:** G1 lag-reduction *below* the PR-3 plateau, joint-vs-per-bond difference measurable.

## 7. Validation plan

All density gates judge **trajectory RMS over [0,2.8]**, never a single time (t=2.6 is a curve-crossing that inverts the D-ladder — `notes/methodology/revival-validation.md:13-22`), and the headline oracle is **`exact_density_finite` (CTM-free)**, never `measure_ctm_pepskit` (χ=8 flatters by ~the signal size). The instrument that drives truncation (the environment) and the instrument that judges (the exact oracle) are deliberately different.

| Gate | What it checks | Quantitative pass/fail | Instrument |
|---|---|---|---|
| **G0** | Frozen-state exact-best-rank-D vs simple update (settles C-vs-D) | (ii)≈(i) within ~ε while (iii) lags → bare Stage C dead; else Stage C has headroom | `exact_density_finite` + exact star gate, 4×4, D2/D4 |
| **G1** | Fixed-D lag reduction | `rms_envaware[D] < rms_SU[D] − margin` for D=2,3,4; D-monotonicity preserved (`rms3≤rms2`, `rms4≤rms3`). Margin > the ~5e-4 SVD-tie null baseline | `exact_density_finite(16)`, extend `test/test_d_convergence.jl:67-122` |
| **G2** | Anti-flattering: monotone in χ (CTM) / patch radius (NTU) | CTM-free trajectory RMS monotone non-increasing AND plateaus (`\|rms(χ_max)−rms(χ_prev)\| < 1e-3`), η≤tol at plateau | `exact_density_finite`, NOT CTM |
| **G3** | Energy-drift tripwire (orthogonal, gauge-invariant) | `\|E(t)−E(0)\|/\|E(0)\| < 5e-3` AND **not** small while density RMS is bad (over-damping signature → FAIL) | **BUILD** `exact_pxp_energy_density_finite(::PXPIPEPSState)` — mirror `_native_dense_state`, apply star H per center, ≤16 sites, ledger-independent (does NOT exist; native energy paths are mean-field or CTM-contaminated) |
| **G4** | ED-free reverse-echo R(t) | `R_envaware(t) ≤ R_SU(t)` across [0,2.8]; isolate truncation from Trotter by varying dt (R→0 as dt→0 ⇒ Trotter; residual ⇒ truncation irreversibility) | **Gauge-invariant** exact-contraction distance (4×4) or frozen-ED endpoint (6×6); replace the mean-field metric in `validate_pxp_reversibility` (`PXPValidation.jl:640-663`); reuse `reverse=true` (`PEPSKitEvolution.jl:154-168`). Caveat: the reverse step is the exact inverse only if `[P,U]=0` AND no truncation (P6) |
| **G5** | Multi-cell anti-overfitting | Lag reduction holds on BOTH 4×4 (native exact + `neel_to_revival_4x4.json`) AND 6×6 (frozen ED `neel_ed_krylov_6.json`, density-only, sample on its exact grid). Compare to ED wherever it fits, never only to SU | exact oracle (4×4); frozen 6×6 ED (D2/D3 only — no native exact oracle >16 sites; D≥4 at 6×6 is UNVERIFIED, treat as such) |
| **G6** | Snake-MPO path-independence (Stage D only) | Two snake orderings WITH inter-tensor truncation agree to tol | required before M7 is used to apply (not store) the gate |
| **G7** | Sublattice Z₂ tripwire (new, cheap) | Sublattice density imbalance stays ≤ bare-SU baseline (exact dynamics conserves it) | `exact_density_finite` per sublattice |

**Self-deception guards.** (a) χ-flattering is now *inside* the evolution for CTM-FU — G1/G2 measure with the CTM-free oracle so a contaminated `N` cannot score itself well. (b) Energy non-conservation masking: a too-aggressive clamp damps the state, looking energy-stable while killing the revival — G3 requires energy-bounded AND density-improved *jointly*. (c) Gauge-invariant-but-wrong: energy conservation and reverse-echo are *necessary not sufficient* (any unitary-ish evolution passes) — the exact-oracle comparison (G1/G5) is the only sufficiency check and must always run where it fits.

## 8. Risks & open research questions

- **[RESEARCH-GRADE] The indefinite-metric well-posedness.** Resolved *for the norm-fidelity objective* (PSD, §3) — but this rests on the code reading of `cost_function_als`/`bondenv_fu`, not on an independent proof for *this* model. **Open:** measure the eigenvalue spectrum of the symmetrized CTM `N` at a few revival times to quantify how much negative weight `positive_approx` actually clamps and whether the constrained PXP sector keeps `N` near-PSD. If the clamped-negative fraction is large, even the demoted CTM-FU ceiling is suspect.
- **Per-bond plateau.** Stage-C-done-right (per-arm, uncoupled cross) is predicted to plateau above the D-ladder; whether the residual is large enough to justify Stage D's cost is the central staging question G0+G1 must quantify.
- **Cost is unmeasured (P4).** Every "in-core at D≤4" / "acceptable cost" statement is back-of-envelope; the dt-vs-staleness coupling to order-2 Trotter cancellation is entirely unmodeled.
- **`E₊` aliasing on 4×4 (P5).** If a CTM plus-star env self-aliases on the only cell with an exact oracle, Stage D's CTM variant cannot be validated there — forcing the NTU finite patch on 4×4.
- **Ephemeral evidence.** The 95–350× error/discard ratios live only in `/tmp/probe{2,4}.jl`; the D4 run was cut at t=2.4. PR-0 must persist them. The M7 rank-2 and path-independence claims lack committed contract tests.
- **No native exact oracle >16 sites.** 6×6 D≥4 has no fast trusted contractor (dense OOMs ~244 GB at D5; `:boundary` is legacy-only, ~46 min/density at D5, OOMs at the revival peak). G5 at 6×6 is D2/D3 only.
- **Open questions:** minimum NTU patch radius that beats SU at fixed D? does joint 5-tensor ALS beat per-arm enough to justify cost at the entangled peak? is the revival-rise lag density-specific or also in energy/fidelity/Loschmidt echo (changes the diagnosis if density-only)? is t≈2.6 on 4×4 the right target or a finite-size artifact (validate via G4/G5 at 5×5/6×6)?

---

Key file:line anchors for the implementer: truncation hook `PEPSKitStarUpdate.jl:207-232` (replace `_split_bond` at `:322`; reuse `Q` at `:285`; deprecated `_bond_trunc` at `:196-202`); gate-in-environment `PEPSKitStarUpdate.jl:280,290,299`; gate projection `PXPModel.jl:200-211`; residual mislabel `PEPSKitMeasurements.jl:493,497-499`; native context (no η gate) `PEPSKitObservables.jl:181,199-210`; exact oracle `PEPSKitObservables.jl:474-553`; trajectory gate `test/test_d_convergence.jl:67-122`; M7 `PEPSKitStarSnake.jl:78-106`; PEPSKit primitives `bondenv_fu` `benv_ctm.jl:24`, `positive_approx` `gaugefix.jl:14`, `calc_convergence` `ctmrg.jl:191-207`.

---

## Addendum — independent post-synthesis verification (2026-06-04)

The synthesis was spot-checked against the installed PEPSKit 0.7.0 source
(`~/.julia/packages/PEPSKit/zz9BS/`) and the working tree. **Confirmed:**

- **The norm-fidelity-is-PSD pivot is correct (§3, P0).** `bondenv_fu`
  (`benv_ctm.jl:24`) builds the bond environment as the bra–ket double layer
  `XX = X'X`, `YY = Y'Y` closed against the CTMRG corners/edges — i.e. the norm
  environment ⟨ψ|ψ⟩, PSD by construction in the χ→∞ / exact-patch limit.
  `positive_approx` (`gaugefix.jl:14`) does `eigh((N+N')/2)`, a global sign flip,
  clamp-negatives-to-zero, sqrt — and PEPSKit's own source comment states the
  negativity arises *"e.g. obtained approximately from CTMRG"*, i.e. it is a
  **finite-χ approximation artifact, not fundamental real-time indefiniteness.**
  This confirms the migration spec's "indefinite metric is ill-posed" framing is
  a category error (it conflates the norm env with the unused energy env), and is
  why the plan picks the NTU finite-patch (PSD-by-construction, no
  positivize-by-fiat) over CTM-FU as the primary engine.
- **All cited PEPSKit primitives exist at the cited locations:** `cost_function_als`
  (`als_solve.jl:175`), `bond_truncate` + `ALSTruncation` (`bond_truncation.jl:58,21`),
  `FullEnvTruncation` (`fullenv_truncation.jl:26`), `calc_convergence`
  (`ctmrg.jl:191,202`) — the real η judge P1 needs, which is absent from
  `leading_boundary`'s returned `info`.

**Still UNVERIFIED (load-bearing — this is precisely why PR-0/G0 exists):**

- The **95–350× error-to-discard-weight ratio** that kills bare Stage C lives only
  in `/tmp/probe2.jl` (D=2) and `/tmp/probe4.jl` (D=4): the scripts survive but
  their stdout was never captured, and the D=4 run was cut at t=2.4. The key
  numbers are recorded in §2 and above, but the prescribed frozen-state
  exact-best-rank-D diagnostic (G0) has **not been run** — PR-0 both runs G0 and
  commits the probe artifacts into the repo.
- One completeness-critic gap the synthesis under-weighted: **log-norm / norm
  bookkeeping under non-unitary real-time env-weighted truncation.** `evolve!`
  accumulates `add_log_norm!` (`PEPSKitBackend.jl:94`) assuming the bare-SVD
  normalization; an env-weighted truncation changes the norm differently. Density
  is a ratio (partly robust), but energy and Loschmidt echo (gates G3/G4) are NOT
  — fold a norm-consistency check into prerequisite P3.