# Fan-out review verdict — CTM-upgrade workflow (2026-06-07)

**Method.** 69-agent fan-out adversarial review of
`docs/superpowers/plans/2026-06-07-ctm-upgrade-workflow.md`: recon ×3 → 6-lens review
(physics / PEPSKit-API / feasibility / sequencing / validation / strategy) → 3-way
refute-by-default verification → 3+1 decision panel → 2 completeness critics. 55 findings
→ 45 material → 18 verified → **10 confirmed blockers/majors**. Every load-bearing fact
below was then **independently re-confirmed against the files** (not taken on the agents'
word).

## Verdict (panel-unanimous, confidence 0.85): do NOT start with PR-0/G0

Make the first slice a combined **de-risking spike** and **un-lock engine decision #1**
(exact-finite-cluster) pending its outcome.

## Verified blockers

1. **The exact-cluster env is NOT "the oracle contraction reused" [F16, blocker].**
   `_native_dense_state` (`PEPSKitObservables.jl:497,504`) is **single-layer** — physical
   legs OPEN. The `{2,2}` bond env is **double-layer** (bra-ket, physical legs traced, the
   two bond-end virtual legs open): a strictly heavier `:boundary`-class object. The
   state-doc itself says this (`...state.md:147-148`); the workflow's "reuses the
   already-validated oracle contraction / zero new infrastructure" (workflow.md:131,151,289)
   is therefore false.
2. **The lighter single-layer oracle already fails at D=4 [F16; re-verified].**
   `g0_frozen_step.out`: D=2 completes all four times; **D=4 killed (signal 15) after only
   t=1.6, 21.3e9 allocations, GC 7798**, inside `exact_density_finite`
   (`PEPSKitObservables.jl:550`). The heavier double-layer env on a full D=4 trajectory — the
   M3 headline proof — is unproven and may be infeasible. (Native `exact_density_finite` is
   `:dense`-only; it rejects the memory-safe `:boundary` fallback — `:537-542`.)
3. **The native D-ladder is FLAT at D≥3 [F1, F2; re-verified].**
   `native_d34_trajectory.out`: native D=3 ≡ D=4 to **6 digits at every time** (identical
   `chunk_trnerr` too); on the rise native lags ED *more* than legacy (`n_native` vs
   `n_leg4`). A1's monotone ladder (D2 2.6e-2 → D3 1.8e-2 → D4 9.8e-3,
   `revival-validation.md:16-22`) is **legacy-backend only** — not the native path PR-3
   repairs. A flat native floor is the **Trotter/scheme signature**, and **no dt-control
   exists** (all probes hardcode dt=0.02) to rule Trotter out. The design-doc itself expects
   "a plateau above the simple-update D-ladder" (`...evolution.md:95`).
4. **The headline gate G1 measures the WRONG backend [F34; re-verified].**
   `test_d_convergence.jl:97` builds `checkerboard_square_ipeps` + `TrotterParams` +
   `evolve!` = **legacy ITensors bare-SVD**, not native `project_star_pepskit!`/
   `evolve_pepskit!`. The native D=2 blowup A6 promises to fix "as a bonus" is **native-only**;
   this test never exercises it. The ED array + thresholds are calibrated to the legacy ladder.
5. **The engine choice INVERTS the plan-of-record [F42; re-verified].**
   design-doc:41 — "**NTU-style PSD finite-cluster environment as the primary engine**";
   :95 PR-3 = NTU; :96 PR-4 = cluster. The workflow locked **exact-cluster** (state-doc:134
   "NEW idea this session … Not yet built", admittedly 4×4-only, ships **no** production
   artifact) as engine #1, while **NTU** — required for 6×6 (G5) — has **zero prototype**
   anywhere (grep: no patch-env code in `src/` or `scripts/`).
6. **Two overstatements [F37, F12; re-verified].** `exact_pxp_energy_density_finite`
   **exists and is exported** (`Observables.jl:476,926`; `SquarePXPDynamics.jl:117,288`) —
   G3/A8 needs a `::PXPIPEPSState` method, not green-field (workflow.md:184 wrong). And A2's
   "PSD ⇒ no `positive_approx` / sidesteps the CTM-trust apparatus" is false: `fixgauge_benv`
   consumes `Z` from `positive_approx`, and cond(Z'Z)=3.4e20 on the rise (`ctm_bond_proto.out`).

## Recommended redirect — the de-risking spike (one slice; de-risks BOTH engines)

- **N4 (composition).** Extend `scripts/g0/ctm_bond_proto.jl` to feed a **hand-built
  exact-cluster `{2,2}` double-layer env** (X/Y isometries absorbed, the two q-reduced bond
  legs open, physical legs traced) through `PEPSKit.bond_truncate` on ONE horizontal bond,
  with a hard `space(hand_env)==space(bondenv_fu(...))` + dual-flag assert and a
  no-op-at-maxdim=BIG self-check. *(The only validated prototype uses the **CTM** `bondenv_fu`
  env the engine demotes, not the exact-cluster env PR-3 needs; its leg-glue is an admitted
  guess; `ctm_star_truncate.jl` stalled before output.)*
- **N2 (feasibility).** Build that same env once at frozen **D=2 and D=4** (t≈2.0–2.4, below
  the t≈2.6 peak), measure wall-clock + peak RSS, extrapolate to the per-step × per-star ×
  per-bond × per-arm-rebuild trajectory cost.
- **Controls.** One D=4 4×4 Néel trajectory [0,2.8] at **dt=0.01 vs 0.02** (Trotter
  attribution) + the **native** D=2/3/4 trajectory-RMS ladder.
- **Decision rule.** Commit exact-cluster only if N4 composes ∧ N2 fits in memory at D=4 ∧ the
  lag is dt-stable ∧ the native ladder is non-flat. Otherwise → **NTU-first** (NTU shares the
  same `{2,2}` composition surface, so N4 de-risks it too).

## Factual corrections to the workflow doc (apply regardless of the spike outcome)

- Strike "reuses the validated oracle contraction / zero new infrastructure / runnable now"
  (workflow.md:131,151,289). State the honest cost class: double-layer `:boundary`
  (~46 min/density, ~244 GB at D=5 entangled peak, `improvement-roadmap.md:27-53`).
- Re-author **G1** as a **native** trajectory test; FIRST prove the native dense oracle
  survives a full D=2 *and* D=4 [0,2.8] run before crediting the gate.
- Add a **dt-convergence gate** + the native D-ladder before crediting A1; restate the
  env-fix success target as beating the **native** floor, not legacy 9.79e-3.
- Re-size **G3/A8** to "add a `::PXPIPEPSState` method to the existing
  `exact_pxp_energy_density_finite`"; add a legacy-vs-native equivalence cross-check.
- Drop A2's "needs no `positive_approx`" selling point; keep `positive_approx→fixgauge_benv`
  on the critical path + a cond(Z'Z) tripwire in PR-3's gate set.
- Close the **cross-observable** question (energy/Loschmidt) before asserting A9
  (density-only adequate) — flagged open at `...evolution.md:124`.
- **Decouple PR-1** (real CTM η) from the PR-3 critical path: P1 is "CTM-FU path only"
  (design-doc:81); the exact-cluster/NTU env bypasses CTMRG. Keep P0 + de-trap as parallel
  low-risk work.

## Open caveats (from the completeness critics — NOT yet resolved)

- "signal 15" is SIGTERM (did-not-complete / thrashing), not provably a strict OOM — **N2
  pins it**.
- **Norm / log-norm bookkeeping** under non-unitary env-weighted truncation
  (`add_log_norm!`, `PEPSKitBackend.jl:94`) is unaddressed by both the workflow and the
  review (design-doc:163-168). Density is a ratio (robust); energy / Loschmidt are not.
- Whether the t≈2.6 4×4 revival is a **finite-size artifact** vs a true scar revival is still
  open (design-doc open-questions).
- Whether an energy/Loschmidt cross-observable is even **well-defined** under non-unitary
  env-weighted steps is unchecked — the A9 fix may not be measurable as posed.

---

## Spike results (2026-06-07, `derisk-spike-exact-cluster` workflow)

The de-risking spike ran (design → spec-review → implement → result-review → decide; 11 agents).
A competing 24-thread `dsfscript.jl` saturated first-call JIT compilation, so the heavy timed runs
were killed mid-compile and several numbers are honestly **UNMEASURED**. File-verified outcomes:

- **N4 composes at the SHAPE level (value-correctness UNVERIFIED).** The spec-review caught a real
  dual-flag bug (`space(X,2)`/`space(Y,4)` are *both* dual; the naive recipe gives
  `codomain≠domain` and would `SpaceMismatch`). The **fixed** hand-built exact-cluster `{2,2}` env
  has space `(ℂ^16⊗ℂ^8)←(ℂ^16⊗ℂ^8)`, flags `[0,0,1,1]`, **byte-identical to the `bondenv_fu`
  ground truth** (`ctm_bond_proto.out`) that already passes `bond_truncate`. **⇒ the SpaceMismatch
  blocker is resolved.** But the no-op-at-BIG / D-truncation / negative-control gates **never
  executed** (cold-compile timeout), so env VALUE-correctness is unproven. `n4Composed:yes` =
  shape-only (reviewers flagged the overclaim).
- **N2 UNMEASURED at every D.** The double-layer env-build wall/RSS — the genuinely-new feasibility
  number — was never reached. Not an OOM, not a fit-confirmation.
- **CORRECTION to blocker #2 above (the "D=4 OOM"): FALSE ALARM.** A single `exact_density_finite`
  at D=4,t=2.0 completes in **16.4 s / 2.5 GB** (754 GB node). `g0_frozen_step`'s `signal 15` was
  allocator thrash (~33 dense calls/time-point, no GC between → 21.3e9 allocs → SIGTERM), not a
  single-state memory wall; `native_d34`'s D=4 density was correct. **The single-layer dense oracle
  is not the D=4 feasibility blocker.** (16.4 s/2.5 GB is a spec-review figure, not re-measured at
  decision time; the sig-15-is-SIGTERM caveat stands.)
- **DECISION-DOMINANT: the native ladder is FLAT.** native D3≡D4 to 6 digits at every time; native
  *leads* ED by ≈0.046 on the rise (the discrepancy is native-above-ED, not a lag). The ladder is
  **already bond-dimension-converged**, so **no environment engine (exact-cluster OR NTU) can move
  it** — both only improve truncation quality at fixed Trotter scheme.

### Engine call (panel-unanimous, conf 0.82): get the dt-control data FIRST

Do **NOT** commit exact-cluster; do **NOT** flip to NTU-first. Both are downstream of the question
the flat ladder forces: **is the native-vs-ED rise discrepancy a Trotter artifact or an env-fixable
deficiency?** The **dt-control** (`scripts/g0/spike_dt_control.jl`: D=3 4×4 Néel, dt=0.01 vs 0.02
over [0,2.8], rise-region drift) decides it:
- **dt-STABLE** (rise drift ≪ 1e-2) ⇒ environment/representation ⇒ an env engine is worth
  committing → then finish N4-content + N2-D4 to choose exact-cluster vs NTU.
- **dt-SENSITIVE** (drift ~ the 1e-2 signal) ⇒ Trotter ⇒ **the PR-3 env thesis must be rethought**
  (the fix is a smaller dt / higher-order Trotter, not an environment) — no truncation engine helps.

**Status:** dt-control launched 2026-06-07 (backgrounded + niced; ~3–4 h, RSS-bounded). N4-content
and N2-D4 held until it returns.

---

## dt-control review (2026-06-07) — ⚠️ ED REFERENCE BUG; the engine thesis is REVIVED

The dt-control ran (exit 0) and its auto-verdict said "dt-SENSITIVE → Trotter." A fan-out review of
that result (4 interpreters → 3 adjudicators → synthesis) **overturned both the auto-verdict and the
spike's "flat ladder → no engine" conclusion** by finding a data bug, which I then **re-verified
computationally against the canonical ED artifact**:

- **CONFIRMED BUG — the `ED02` reference is the true ED shifted +0.2 on the rise.** In *both*
  `scripts/g0/spike_dt_control.jl:17-18` and `scripts/g0/native_d34_trajectory.jl:11-12`,
  `ED02[t] == canonical[t−0.2]` for **every** t≥1.6 (verified exact vs `artifacts/neel_to_revival_4x4.json`).
  The collapse branch (t≤1.4) was spliced from a D=4 iPEPS output, not ED.
- **CONSEQUENCE — the "native LEADS ED by ≈0.046 / flat ladder ⇒ no engine" claim above is RETRACTED.**
  Against the **true** ED, native dt=0.02 D=3 RMS = **8.63e-3** and it **LAGS** ED by **−0.0094** (the
  documented small representation lag) — not a +0.046 lead. The +0.046 "lead" was 100% the shift.
- **The dt-anomaly is real and anti-Trotter.** The order-2 palindrome + `_nsteps` are verified correct,
  so "smaller dt → closer to ED" is the right expectation — yet dt=0.01 D=3 RMS = **8.10e-2** (~9.4×
  *worse*; revival destroyed, 0.349 vs ED 0.483). Even in the collapse region dt=0.01 is 1.6–3.1× worse
  (order-2 demands ~4× better) and Richardson **diverges** ⇒ both dt are **outside the asymptotic
  regime**; a competing error *grows* as dt→0.
- **That competing error is truncation-in-a-bad-environment, not bond budget.** Same mechanism as the
  rel_floor-independent native D=2 blowup at the identical t=1.6→1.7; `err/truncerr` decoupling 95–1156×.
  The D3≡D4 "flat ladder" is likely a `rel_floor=1e-3` artifact (the floor cuts beyond D=3) + bad-subspace
  selection, **not** true bond-convergence. dt=0.02's good agreement is partly a fortuitous
  Trotter↔truncation-gauge cancellation; dt=0.01 breaks it.

**Net (panel conf 0.85, my numbers re-verified):** the engine-killing premise is dead, so the
**environment-engine thesis is REVIVED** — the env-is-the-bottleneck evidence (truncerr/lag decoupling,
rel_floor-independent D=2 blowup, 26× CTM-aware-vs-bare-SVD) is **independent of the bug and survives**.
BUT engine commitment stays **DEFERRED**: there is no dt-converged baseline, and the Trotter knob
currently dominates (dt 0.01-vs-0.02 differ by up to 0.13 on the rise, ~13× the ~1e-2 signal an engine
targets). The D-ladder lag numbers (D2 2.6e-2 … D4 9.8e-3) are all single-dt=0.02 and partly fortuitous;
**9.8e-3 cannot be quoted as the engine target until reproduced at dt-convergence.**

**Next experiment:** a **dt-converged (near-lossless) baseline vs the TRUE ED** — dt-ladder at rel_floor=0.
Feasibility caveat I flagged (not in the panel's spec): truly lossless likely needs D≥5, where the dense
oracle is heavy (~244 GB at D=5 per `revival-validation.md`); D=4 dense is cheap (16 s / 2.5 GB).
**To fix:** replace the hardcoded `ED02` in both scripts with the canonical `ed.density` loaded from
`artifacts/neel_to_revival_4x4.json`. **DONE 2026-06-07** (both scripts now carry the true ED).

---

## Scheme audit (2026-06-07) — NO bug; the instability IS the environment ✓ (engine thesis founded)

Read-only fan-out audit (6 suspect paths → refute-by-default verification → decide) of whether the
dt=0.01 anti-Trotter blowup is a fixable scheme bug. **Verdict (no candidate bug survived): NO scheme
bug — all 6 non-truncation paths are correct-as-is; the blowup is INHERENT bad-subspace truncation
against the severed-loop λ² mean-field environment** (a fixed-point step-count feedback). Empirically
ruled out (agents ran isolated checks):
- **Normalization** — `normalize!(T,Inf)` factors out as a single global scalar (bit-identical bond-SV
  ratios under injected scalars; 5.6e-17 toggle agreement) and is *necessary* (disabling overflows at
  t=0.20); the dense oracle is a normalization-invariant Rayleigh quotient (value/normsq). Ruled out.
- **√λ weight floor** — never fires (smallest weight ~4.3e-3, ~10 orders above the ~1.8e-12 threshold);
  absorb/deabsorb exact round-trip on disjoint slots. Ruled out.
- **Gate P·U** — exactly-zero forbidden leak ([P,U]=0). **_split_bond / convention-B gauge / oracle** —
  correct-as-is.
- **Discriminators all point at the environment, not a bug:** the blowup is anti-convergent in dt
  (dt=0.01 ~9.4× *worse* vs Trotter's expected ~4× *better*), D-independent (D3≡D4 to 6 digits), and
  rel_floor-independent (0 vs 1e-3 identical). A bug would be sensitive to ≥1 of {disable-toggle, D,
  rel_floor}; none moves it.

**Consequence — the env-engine thesis is now WELL-FOUNDED** (conclusively the environment, not a bug),
**and the dt=0.01 D=3 blowup is the cleanest, largest test case for the engine**: the revival is fully
destroyed (peak 0.349 vs ED 0.483), a big unambiguous signal — far better than the original small ~1e-2
lag. The decisive engine test is now a clean go/no-go: **swap ONLY the truncation environment** (bare λ²
→ loop-aware exact-cluster/NTU — which the spike showed *composes*) at the SAME dt=0.01, D=3; if the
blowup vanishes with every audited path byte-unchanged, the environment is conclusively the cause AND
the engine works. **Do not chase dt→0** (it makes the revival worse). Optional cosmetic-only: fold the
`normalize!(Inf)` scale into `add_log_norm!` (does not change any observable).

---

## PR-3 env prerequisites (2026-06-07) — exact-cluster env WORKS but is trajectory-INFEASIBLE

Ran the spike's env scripts to completion (no timeout, in the main loop). Two results:

- **VALUE-CORRECT ✓** (n4, D=2 frozen, one horizontal bond, `spike_n4_content.out`): no-op-at-BIG dev
  **6.5e-10**; D-truncation fid **1.0**, |dn| **4.30e-9**; and **~92× more accurate than bare-SVD**
  (|dn_bare| 3.97e-7) on the bond — the loop-carrying env demonstrably helps. cond(Z'Z)=9e20 →
  `fixgauge_benv` needed. Negative control (global-conj env) → NaN (env matters; weak control). The
  50% raw-tensor diff vs the χ=16 CTM env is expected (exact ≠ χ=16-approx; observables validate it).
- **TRAJECTORY-INFEASIBLE ✗** (n2, `exact_cluster_env_cost_d2.out`): one D=2 double-layer env build =
  **127 s / 8.1 GB** (12× the single-layer oracle). Per-step = 64 builds → full trajectory
  **≈ 316 h ≈ 13 DAYS at D=2** (worse at D=3; the builder is also D=2-hardcoded — `@assert DX dim 16`
  errors at D=3 where q-legs are 24).

**Consequence — the exact-cluster env is the gold-standard accuracy REFERENCE, not the production
trajectory engine** (exactly the design-doc's NTU-primary call, now empirically vindicated). The
decisive dt=0.01 D=3 trajectory test **requires the cheap NTU finite-patch {2,2} env** (local
~9–13-site double-layer plaquette-ring, not a full-torus contraction — should be ≲ seconds/build),
which is **unbuilt**.

### Path forward
1. Build the **NTU finite-patch {2,2} env** (cheap, local; the design-doc's primary engine).
2. Validate NTU against the exact-cluster gold standard on a few frozen bonds (NTU → exact as the
   patch grows; reuse `build_exact_cluster_bondenv` as the reference).
3. Run the decisive **dt=0.01 D=3 env-aware (NTU) trajectory** vs bare — does the blowup vanish?

The exact-cluster work was not wasted: it proved a loop-carrying env is value-correct + 92× more
accurate, and it is the reference to validate NTU against.

---

## CORRECTION + reframe (2026-06-07) — exact-cluster IS feasible (cold-compile artifact); test at D=2

**RETRACTION (verified):** the "127 s / ~13 days infeasible" premise was a **cold-compile artifact.**
`exact_cluster_env_cost_d2.out`: the FIRST D=2 build (t=2.0) = **127.083 s**, the SECOND (t=2.4, same
D=2) = **2.755 s** — a 46× drop, ~98% first-call ncon JIT. **Warm exact-cluster build ≈ 2.8 s → full
D=2 trajectory ≈ 6.9 h (FEASIBLE).** My earlier "infeasible" conclusion (extrapolated from the cold
first build before the warm number landed) is wrong; the pivot to NTU was made on a bad number.

**NTU P6 (built this round, `reviewer_ntu_probe.out`):** the maximal-patch==exact equality gate
**passes** (builder correct); P6-trace warm build 1.338 s vs warm-exact 7.04 s (5.3× faster, not 100×);
better-conditioned (cond 4.3e13 vs 9e20). BUT the single-frozen-bond |dn| ladder is **near-lossless
noise** — bare 6.1e-6, exact 1.2e-6, NTU **6.9e-8** (NTU "beats" the *exact* env, which is unphysical),
env rel-dev P6-vs-exact = **0.949** (nearly orthogonal). ⇒ the single-bond test does NOT discriminate
env quality (the truncation barely matters at one frozen bond), so NTU's quality is **unvalidated**.
Consistent with the audit: the blowup is **cumulative — only a trajectory discriminates.**

**Reframed path (cleaner than NTU-at-D=3):** run the decisive trajectory test at **D=2 with the
validated exact-cluster env** (~7 h), where native D=2 bare **blows up** (0.143→0.396 at t=1.6→1.7,
`native_relfloor_test.out`). The env-aware-vs-bare integration is **env-agnostic** (every env is a
{2,2} `BondEnv` fed to `bond_truncate`), so build it once behind a flag (`:bare | :exact_cluster |
:ntu`), run `:exact_cluster` at D=2 → does the blowup vanish? This settles "does the engine work" with
**no env-quality uncertainty and no D=3 generalization.** NTU becomes the later cost optimization once
the engine is proven (and the D≥3 / >4×4 production engine).

---

## CORRECTION² (2026-06-07) — integration VALIDATED; premise re-grounded against the TRUE ED

The env-agnostic 4-bond integration workflow (`wz3wur5p5`, 11 agents) returned
**`integrationValidated: true` (conf 0.86)** AND a Decide-phase **premise blocker**. Both verified
against raw data here. Net: the integration is sound; the *pass/fail criterion* was phantom-grounded
and is now re-stated against the canonical ED. The decisive test is **real, not a phantom** — the
Decide agent itself overcorrected.

### Integration: VALIDATED (all 4 star bonds incl. both vertical)
`scripts/g0/star_env.jl` (`build_star_bond_env`, `project_star_env!`) — D-generalized, 4-direction
{2,2} env via rotate→`_qr_bond`→true-geometry-id emission. Gates (frozen native D=2, t=1.8):
- **G1 vertical env-VALUE oracle** (transpose-oracle, the discriminating gate; no-op-at-BIG is
  env-blind plumbing and provably toothless for a leg transpose): up=5.99e-14, down=2.08e-14. PASS.
- **G2 per-dir single-bond no-op-at-BIG:** R=5.6e-15 U=3.4e-13 L=2.7e-14 D=1.9e-15, all <1e-9. PASS.
- **G3 full-star no-op-at-BIG:** ALS(50)=4.97e-9 (>1e-9, ALS iterative floor ~1e-12/arm × 4);
  **`trunc_alg=:full` (FullEnvTruncation)=1.05e-10 <1e-9 PASS.** ⇒ trajectory engine MUST be `:full`,
  not the default `:als` (ALS leaves a ~5e-9/step no-op floor that contaminates a science result).
- **G4 5-step D=2 `:exact_cluster` trajectory:** all finite, density differs from bare 2.9e-4→2.9e-3.
- **`:ntu` vertical is NOT gated** (hand-rolled patch, no NTU==exact for up/down). Use `:torus` /
  `:exact_cluster` for any trajectory; defer the `:ntu` leg until its vertical patch is gated.

### Premise re-grounding (verified against raw .out + the canonical artifact)
1. **Phantom ED — CONFIRMED.** `native_relfloor_test.jl:12` `ED02` was the canonical ED **shifted
   +0.2** (`ED02[1.6]=0.10056 == canon[1.4]`). FIXED this round to the true values
   (`1.6=>0.15558, 1.8=>0.24370, 2.0=>0.33735, …`).
2. **But a genuine D=2 blowup vs the TRUE ED — ALSO CONFIRMED** (`native_relfloor_test.out`,
   rel_floor=0, **min_bd=2 throughout — NOT a D=1 drop**): D2 tracks ED to t=1.6 (0.145 vs ED 0.156),
   then **jumps to 0.403 (t=1.7), 0.428 (t=1.8)** — true ED only 0.244, a **+0.18 overshoot** — stays
   elevated through t=1.9 (0.410), then crashes below ED by t=2.2 (0.307 vs 0.416). A real, sustained
   ~0.18 excursion. The Decide agent **understated** this as a "single-point rel_floor spike."
3. **D=2→D=3 REMOVES the blowup; D=3≡D=4** (`native_d34_trajectory.out`, native column is
   label-agnostic physics): D=3 and D=4 identical to 6 digits — smooth, no spike, t=1.8=0.2336 vs ED
   0.2437 (~1% undershoot). So bond-dim **2→3 matters** (kills the spike); **3→4 is converged**. The
   Decide agent's "bond-dim does NOT change the rise" **conflated D=3≡D=4 with D=2≡D=3 — overstated.**

### Net — the env thesis stands, sharpened into two ED-grounded effects
- **(A) D=2 blowup** (bond-sensitive 2→3): the PRIMARY decisive target. Does `:exact_cluster` at D=2
  suppress the +0.18 t=1.7–1.9 overshoot and pull D=2 onto the D=3/D=4 (ED-tracking) curve? I.e. does a
  loop-carrying environment let a SMALL bond dim behave like a larger one — proving *environment, not
  bond-dim* drives the D=2 instability. Clean, falsifiable, decisive.
- **(B) residual ~1% bond-CONVERGED lag** (D=3≡D=4 undershoot vs ED): the harder mean-field-ceiling
  question — can a loop-carrying env beat the simple-update environment that bond-dim cannot move?
  Secondary; only meaningful after (A).

### Path forward (re-stated; preconditions from the Decide phase, all now satisfiable)
1. Phantom ED in `native_relfloor_test.jl` — **FIXED** this round.
2. **Cheap onset-intervention probe FIRST** (before any 6–12 h run): take bare D=2 to t=1.6 (just
   before onset), then run bare vs `:exact_cluster, trunc_alg=:full` through ~t=2.0 on the 0.1 grid,
   vs true ED. Isolates the truncation-ENVIRONMENT effect at the disease onset (cleaner control than
   from-t=0, which confounds accumulated history) AND exercises fixgauge on the **real rise**
   (cond~1e20 — the regime gap every frozen-t=1.8 gate left unproven). ~minutes of env builds.
3. **If the probe suppresses the onset jump** → write the full-sweep driver
   (`scripts/g0/decisive_env_trajectory.jl`: sweep all reps in the :serial order-2 palindrome calling
   `project_star_env!(env=:exact_cluster, trunc_alg=:full, patch=:torus)`) and launch the from-t=0
   D=2 trajectory in background (~6–12 h) for the publication curve. Existing scripts call
   `project_star_env!` on ONE center/step — no full-sweep driver exists yet.

---

## DECISIVE RESULT (2026-06-07) — env-aware D=2 SUPPRESSES the blowup, lands on the D=3/D=4 curve

Onset-intervention probe (`scripts/g0/onset_intervention_probe.jl`,
`onset_intervention_probe.out`): bare D=2 evolved to the disease onset t=1.6 (rel_floor=0,
min_bd=2), then forked — same :serial order-2 schedule, swapping ONLY the truncation engine.

| t   | bare (native) | env :exact_cluster :full | D=3≡D=4 | true ED |
|-----|---------------|--------------------------|---------|---------|
| 1.6 | 0.14480       | 0.14480                  | 0.14615 | 0.15558 |
| 1.7 | **0.40269**   | **0.18052**              | —       | —       |
| 1.8 | **0.42815**   | **0.22321**              | 0.23363 | 0.24370 |
| 1.9 | 0.40987       | 0.26977                  | —       | —       |
| 2.0 | 0.38292       | 0.31652                  | 0.32746 | 0.33735 |

**The loop-carrying {2,2} environment at D=2 removes the +0.18 blowup and reproduces the
bond-converged (D=3≡D=4) curve to ~0.01.** At t=1.8: env 0.223 vs D34 0.234 (Δ0.010) vs ED
0.244; bare 0.428 (Δ vs ED 0.184). This confirms the thesis: the D=2 instability was a
truncation-ENVIRONMENT artifact; an exact-cluster environment lets a SMALL bond dim behave like
a larger one. The residual ~1% env-D2-vs-ED undershoot == the bond-converged mean-field lag
(effect B), a separate/harder problem; effect A (the blowup) is DECISIVELY fixed.

Diagnostics (all 4 chunks): min_bd=2 (genuinely D=2), nfb=0 (fixgauge_benv NEVER fell back —
robust even at cond(Z'Z)→Inf on the rise; CLOSES the "regime gap" the consensus reviewers
flagged — every prior gate ran at frozen t=1.8 cond~1e6), finite=true throughout. bare branch
reproduces native_relfloor_test.out's spike to 4 digits (0.4027/0.4281) => faithful control.
FET "cancel 50" warnings = FullEnvTruncation hitting maxiter=50 at fid≈1−3e-7 (converged-enough
for D=2). Cost: ~30 min/0.1-chunk (640 exact-cluster env builds × ~2.8 s) => full [0,2.8] ≈ 14 h.

**Integration bug found+fixed en route** (only a MULTI-CENTER sweep exposes it; all G1–G4 gates
were single-star, structurally blind to it): `env_truncate_bond!` wrote the truncated tensors to
`peps.A` but NOT the SUWeight ledger. The lossless star grows the bond weight to full size; env
truncation shrank peps.A but left `weights.data[slot]` stale → the next center's `absorb_weight`
threw `SpaceMismatch(ℂ^2 ≠ ℂ^8)`. Fix mirrors project_star_pepskit!:367-368 (store non-dual S at
`_bond_slot`), in BOTH the env and fallback branches of env_truncate_bond!.

### PENDING adversarial control (causal attribution)
project_star_env! grows all bonds lossless THEN re-truncates — so the fix could be the loop
ENVIRONMENT *or* merely the lossless-then-retruncate procedure. CONTROL: add env=:bare
(lossless → bare-SVD truncate, no env) and check whether it spikes like native-bare. :bare
spikes + :exact_cluster smooth => environment decisively the cause. (Cheap: no env builds.)

---

## REVIEW of the decisive result (2026-06-07, workflow wedw4yosa, 10 agents) — claim (c) RETRACTED

Adversarial fan-out review of the onset result + causal control. Verdict by claim:

- **(a) spike = inline-ORDERING — QUALIFIED.** Data holds (both bare variants de-spike ⇒ the
  environment is NOT what removes the spike). But the LABEL "inline-ordering artifact" is
  underdetermined: env=:bare changes TWO axes vs native-bare (different SVD partition — a 2-body
  _combine_ab+SVD vs the sequential 5-body θ split — AND deferred truncation), so smoothness could
  be the partition, not the deferral. Also D=3/D=4 truncate inline yet don't spike ⇒ "inline active"
  isn't sufficient; the real variable is D=2 subspace selection under a bad mean-field environment.
- **(b) bare-SVD undershoots + diverges — SUPPORTED.** Deterministic, monotone undershoot
  (−0.011→−0.074→−0.142 vs ED), lossless round-trip (Δ~1e-16), normalization inert (Δ<5e-8),
  scale-invariant observable rules out decay. Genuine subspace-bias inaccuracy.
- **(c) loop-env "tracks D34 ⇒ accurate production engine" — RETRACTED (refuted as written).**
  D34 is NOT bond-converged (the D3≡D4 "flat ladder" is itself a rel_floor=1e-3 artifact, not
  convergence) and dt=0.02 is NOT Trotter-converged; exact_cluster's error vs **true ED GROWS**
  (−0.0108→−0.0205→−0.0208), it does not lock on. "exact_cluster≈D34" only shows two same-dt,
  same-rel_floor, non-converged curves share a ~0.01–0.02 error floor over a 0.4-in-t window — it
  does NOT validate against physics.

**Corrected claim (the defensible one): EFFECT A only.** The loop-carrying environment removes the
catastrophic +0.18 D=2 inline blowup and lands ~7× closer to true ED than bare-SVD re-truncation at
t=2.0 (env Δ_ED=−0.021 vs bare-SVD Δ_ED=−0.142). Engine commitment as "accurate" stays **DEFERRED**
pending a dt-converged baseline. Do NOT write "validated as the accurate production engine" or
"reproduces bond-converged dynamics." (Supersedes the strong wording in the DECISIVE RESULT section.)

**The decisive test for the strong claim is a dt-CONVERGENCE point, not a longer t-window:** from the
t=1.6 fork, env=:exact_cluster:full at dt=0.01 (ideally +0.005) to t=1.8. If exact_cluster's error vs
**true ED(1.8)=0.24370** SHRINKS as dt→0 (the dt=0.01 point moves up from 0.2232 toward 0.2437) ⇒ the
env converges to true physics, (c) rescued. If it instead tracks D4 and both drift from ED ⇒
shared-dt Trotter↔truncation cancellation, (c) stays refuted.

**Confirmed must-fixes:** (1) DATA HYGIENE — the ED *column* in native_d34_trajectory.out /
native_relfloor_test.out is the stale −0.2-shifted label (scripts already corrected; the .out files
are pre-fix). Regenerate native D3/D4 on the 0.1 grid to t=2.8, ED labeled same-t. (2) cond(Z'Z)=Inf
reading CONFIRMED correct (fallback fires only on a caught exception; cond is diagnostic-only) — but
the peak t>2.0 is unexercised; pre-flight one env chunk near t=2.4 before the long run.

**Path forward (confidence 0.86):** engine = env=:exact_cluster trunc_alg=:full (NOT NTU — vertical
patch ungated). Two runs, in parallel (node idle, 754 GB): (i) **dt-convergence** env fork t=1.6→1.8
at dt=0.01 — decisive for (c); (ii) **effect-A-through-revival** env fork t=1.6→2.8 at dt=0.02
(~6.3 h, 12 chunks) — does env-D2 stay de-blown-up through the revival peak (t=2.6) or develop its own
instability past t=2.0? [0,1.6] carries no env-quality info (native bare already tracks ED there), so
fork-at-1.6 replaces the 14.6 h from-t=0 run.

---

## dt-CONVERGENCE RESULT (2026-06-08, `env_dtconv_dt01_to1.8.out`) — gap is METHOD error; (c) strong stays refuted

The decisive dt-convergence point ran (exit 0). Both forks start from the **identical** bare state
(t=1.6, n=0.144802, baredt fixed at 0.02), so the env evolution dt over 1.6→1.8 is the only variable:

| dt    | env(t=1.8) | true ED(1.8) | gap to ED |
|-------|------------|--------------|-----------|
| 0.02  | 0.223205   | 0.24370      | −0.020495 |
| 0.01  | 0.223225   | 0.24370      | −0.020475 |

Halving dt moves env(t=1.8) by **+2.0e-5** — essentially zero (raw 6-digit values; the earlier 0.22321
was the rounded EXC-dict figure). Order-2 Richardson (n(dt)=n₀+c·dt²): c≈0.067, **dt→0 extrapolant ≈
0.22323**, still **−0.0205 below ED**. So the ~0.021 undershoot is
**method/truncation (bond-dimension) error, NOT Trotter** — it does not close as dt→0. Applying the
review's verbatim rule: the dt=0.01 point did **not** move up toward 0.24370 ⇒ **claim (c) strong
reading stays REFUTED, now with a dt-convergence margin** (not just "tracks D4"; it's pinned).

**Two things are simultaneously true and both matter:**
- **(+) The env is dt-STABLE** — 0.02 and 0.01 agree to 1e-5. This is the *opposite* of native bare
  D=3, which was anti-Trotter (dt=0.01 ~9.4× worse, Richardson diverged). The loop environment removes
  the dt-sensitivity pathology too, not just the blowup. The env evolution is well-posed in dt.
- **(−) The converged limit is ~2% below ED** — the env converges to *its own* D=2-truncated dynamics,
  ~0.011 below the D34 baseline and ~0.021 below true ED. To close that you need **more bond dimension
  (D>2), not smaller dt.** Effect A (de-blowup) is dt-robust; faithful-physics is D-bounded.

**Engine implication:** env=:exact_cluster:full is a *stable, dt-converged* D=2 engine that lands ~2%
below ED — usable as a stable integrator, but not a faithful-physics claim at D=2. The accuracy lever
is now unambiguously **D**, not dt.

### Review of the dt-convergence result (2026-06-08, workflow wbmdbcthf, 5 agents) — UNANIMOUS, conclusion HOLDS

Independent re-derive + 3 refute-by-default angles + forward. All 5 confirm: gap = **method/D=2-truncation
error, not Trotter**; claim (c) strong (fixed-D=2 → true physics) **stays refuted**. No refute angle was fatal:
- **Asymptotic/underpowered (refute-0):** a 2-point dt test can't *prove* asymptotic membership, BUT any
  analytic Trotter fit (order-2 OR order-1) lands at 0.2232, **~3000× short** of the +0.0205 climb needed
  to reach ED. The only escape is a super-singular sub-dt=0.01 feature pushing *up* — yet the one documented
  small-dt pathology (native-D3 severed-loop feedback) pushes the **wrong sign** and is removed by
  construction here (loop-aware {2,2} env replaces the severed-loop λ² env). dt=0.005 would be
  confirmation-only.
- **FET-floor confound (refute-1):** FET hits `cancel 50` (maxiter) not the tol break, but it is stationary
  (Δfid ~1e-13–1e-15) and the floor is **rank-2 representability (dt-independent)** — so this *strengthens*
  the method-error reading. (Caught a stale comment: `star_env.jl:400-401` calls FullEnvTruncation "exact
  (non-iterative)" — it is iterative, maxiter=50. Doc-only.)
- **Patch-vs-D rescue (refute-2):** the strongest possible rescue ("residual is patch error, a bigger env
  reaches ED") is **closed on the code** — `build_star_bond_env`'s `:exact_cluster` branch
  (`star_env.jl`, `for rr in 1:Nr, ccol in 1:Nc`) contracts **all 16 torus cells with periodic ids → it IS
  the maximal exact env for the 4×4 cell.** A perfect env at D=2 still undershoots ⇒ the residual is
  **D-bounded, not patch-bounded.** Clean decomposition (both legs are D, not patch/dt): env-D2(0.223) →
  native-D3(0.234) ≈ 0.010 (pure bond capacity, env already exact) → ED(0.244) ≈ 0.010 (along the D-ladder).

**Net:** the env is a stable, dt-converged D=2 integrator whose converged limit is D-bounded ~2% below ED.
The honest framing is **effect-A + dt-robustness** (NOT "validated production engine" / "reproduces
bond-converged dynamics"). The weaker effect-B claim ("env beats simple-update at fixed D") is the open
question a **D=3 env run** settles.

**Number corrections applied** (verdict unchanged): dt=0.02 = 0.223205 (not 0.22321); dt step +2.0e-5;
c≈0.067; dt→0 extrapolant 0.22323.

### Forward plan (review-unanimous + my checks)
1. **dt=0.005: SKIP** (`should_run_dt005=false`) — lowest-information run; order-2 predicts 0.22323 ±1e-5,
   indistinguishable, for ~2 h. The analytic 3000×-short argument already makes it airtight.
2. **Effect-A window: let it finish** (running, healthy through t=2.0, ~3 h left) — free stability data
   through the t=2.6 revival peak; settles whether env-D2 develops its own instability past t=2.0.
3. **D=3 env fork = the real accuracy test** (`should_run_D3_env=true`) — does the env *climb* toward ED at
   D=3 (effect B alive) or saturate at the D=2-truncated value (effect-A-only)? GATE it behind cheap
   frozen-bond D=3 checks (no-op-at-BIG + truncate-at-D fidelity at q-legs 24/12) BEFORE the multi-hour
   trajectory. Builder D=3-readiness to be verified next (the agents disagreed on where the D=2 q-leg
   hardcode lives — `star_env.jl` vs `spike_n4_exact_cluster_bond.jl`).
4. **Caveat to resolve (cheap_confirmations):** native_d34_0p1grid baseline used `rel_floor=1e-3` →
   the native-D3=0.2336 anchor may be a floor artifact, not true D=3. Confirm before treating 0.2336 as the
   D=3 target for the env comparison.

### D=3 builder gate (2026-06-08, `envaware_integration_d3.out`) — builder VALUE-correct at D=3 ✓

Parameterized `envaware_integration.jl` for D (default 2, backward-compatible) + scaled BIG→256, ran D=3.
The trajectory engine `build_star_bond_env` (star_env.jl) needed **NO edit** — it was already
D-generalized (the D=2 q-leg hardcode lives only in the unused standalone `spike_n4_exact_cluster_bond.jl`).
Frozen native-D3 at t=1.8 = **0.23362680** (matches native_d34's 0.2336 → the D3 anchor is robust). Gates:
- **G0 ✓** rotate roundtrip 0.0, dual flags all true, **qdim=6** at D=3 (not the "24" the PR-3 note guessed).
- **G1 ✓** vertical env-VALUE oracle (transpose-trick, the discriminating gate): up **1.3e-11**, down
  **4.4e-11** (<1e-10). The independently-constructed D=3 env matches the rotate-route env to ~1e-11.
- **G2 ✓** per-dir no-op-at-BIG: **all ~1e-16** (machine precision), cond(Z'Z) 5e4–1.8e6 (well-conditioned).
- **G3 = false (6.76e-9)** but **EXPECTED & benign** — G3 uses the *default* `:als` (the known ALS floor;
  D=2 was 4.97e-9 with :als, 1.05e-10 with :full). The science engine uses `:full`. NOT a D=3 bug.
- **G4 ✓** 5-step D=3 env trajectory finite, differs from bare ~1–4.5e-6 (monotone). Env-D3 runs stably.

**Two anomalies, both resolved benign:** (1) GMRES "stopped without converging" warnings (residual ≤5.7e-10)
are FullEnvTruncation's **inner linsolve** (`fullenv_truncation.jl:245/255`) reaching ~1e-10 — good enough,
flagged vs a stricter tol; at D=3 bond sizes (~36×36) the 2902 ops are ms, no slowdown; G1/G2 prove the
machinery is value-correct anyway. (2) Fixed the confirmed-wrong comment (`star_env.jl:398`) that called FET
"exact (non-iterative)" — both engines are iterative; FET just converges to a tighter floor.

**⇒ Decisive D=3 control launched** (`b22mr6e78`, `env_d3_control.jl`, ~4h to t=2.0): native-bare-D3 vs
env-D3(:exact_cluster:full) from the t=1.6 fork (both rel_floor=0, dt=0.02), vs ED. Decision at t=1.8 —
env-D3 climbs above env-D2's 0.223 toward ED 0.244 (D-lever works) / beats native-D3 0.234 (effect B alive)
/ saturates (effect-A-only, deeper mean-field ceiling).

---

## EFFECT-A WINDOW (2026-06-08, `env_window_dt02_to2.8.out`) — env-D2 survives the FULL revival ✓

env=:exact_cluster:full, D=2, fork@t=1.6, dt=0.02 → t=2.8 (12 chunks, ~29 min each, 7.5 GB peak):

| t   | env-D2  | true ED | Δ(env−ED) |       | t   | env-D2  | true ED | Δ(env−ED) |
|-----|---------|---------|-----------|-------|-----|---------|---------|-----------|
| 1.6 | 0.14480 | 0.15558 | −0.011    |       | 2.4 | 0.46470 | 0.46696 | −0.002    |
| 1.8 | 0.22321 | 0.24370 | −0.021    |       | 2.6 | **0.49275** | **0.48306** | **+0.010** |
| 2.0 | 0.31652 | 0.33735 | −0.021    |       | 2.7 | 0.49243 | —       | —         |
| 2.2 | 0.40227 | 0.41600 | −0.014    |       | 2.8 | 0.48253 | 0.46230 | +0.020    |

**Effect A holds through the ENTIRE revival.** minkept=2, nfb=0, finite at every step to t=2.8;
**maxcond=Inf throughout yet fixgauge_benv NEVER fell back** — this CLOSES the "peak regime (t>2.0)
unexercised" gap every prior gate left open. Where bare-D2 blew up to 0.43 by t=1.7 and was destroyed by
t=2.2 (never reaching the revival), **env-D2 is a stable integrator that reproduces the revival peak at the
right time** (env peaks t≈2.6–2.7, ED peaks t=2.6).

**But D=2 is NOT faithful physics — the error is a ±2% oscillating band, not a uniform undershoot.** env-D2
undershoots on the rise (−0.021), crosses zero near t=2.4, then **overshoots** the peak (+0.010 at 2.6, +0.020
at 2.8). At the peak env-D2 (0.493) is slightly *worse* than rel_floor=1e-3 native-D3 (0.477) — they bracket
ED. Fully consistent with the dt-converged D-bounded method error: a stable, ~±2%-accurate D=2 integrator,
not a physics engine. **This is the headline effect-A result: a loop env turns catastrophically-unstable
bare-D2 into a stable revival-reproducing integrator at ~±2%.**

### D=3 control fork observation — the native-D3 anchor was a rel_floor artifact (CONFIRMED)
First line of `env_d3_control_to2.0.out`: **bare-D3 (rel_floor=0) at t=1.6 = 0.155540 vs ED 0.15558 (Δ −4e-5,
essentially EXACT).** The old "native-D3 = 0.146 (t=1.6) / 0.234 (t=1.8)" anchor was the rel_floor=1e-3 run —
a floor artifact, exactly as cheap_confirmations flagged. True rel_floor=0 bare-D3 is far more accurate, so
the clean effect-B test is **env-D3 vs bare-D3, BOTH rel_floor=0** (the control provides both columns). This
also reframes effect B: at D=3 bare is *already* accurate and stable (no blowup to fix), so the question
becomes whether env-D3 *matches* or *beats* an already-good bare-D3 — not whether it rescues a blowup.

## ⛔ EXACT-CLUSTER TRAJECTORY IS INFEASIBLE AT D=3 (2026-06-09) — killed after 32 h / 70 GB on ONE chunk

The D=3 control **never completed even chunk 1** (1.6→1.7): 32 h wall, RSS exploded 22→**70 GB**, zero
chunks emitted, killed (exit 144). The exact-cluster env is a **full-torus double-layer contraction** whose
cost scales ~ (fat-bond)^width with fat-bond = D². D=2→D=3 takes fat-bond 4→9; worse, **`project_star_env!`
grows every star bond LOSSLESSLY (maxdim=10⁶) before each env build**, so the trajectory builds the torus env
on *grown* bonds (≫3), not the frozen D=3 bonds the gate used. The frozen-bond D=3 gate ran at ~50 s/build
(feasible); the grown-bond trajectory build is hours/build and GB-scale → **infeasible.**

**This empirically vindicates the design-doc's architecture call.** The exact-cluster env was always the
**4×4-D=2 gold-standard accuracy REFERENCE**, never the production trajectory engine (design-doc:41 makes the
NTU local-patch the primary engine for exactly this reason; the full-torus contraction does not scale in D or
system size). Effect A (the headline) is fully established at D=2 with this reference. **Effect B and any
D≥3 / >4×4 dynamics REQUIRE the cheap LOCAL env** (NTU finite-patch {2,2}: a ~9–13-site double-layer
plaquette-ring, no full-torus contraction) — the design-doc's primary engine, which is built (`ntu_bondenv.jl`)
and passes the maximal-patch==exact gate, but whose **vertical (:up/:down) patch is ungated** and whose
single-frozen-bond quality is unvalidated (only a trajectory discriminates — see the CORRECTION² /
PR-3-prereq sections above).

### Strategic fork (the CTM-upgrade goal re-pointed at the scalable engine)
The exact-cluster work delivered its purpose: it PROVED a loop-carrying env (a) is value-correct (G1 1e-11,
G2 1e-16, at D=2 *and* D=3), (b) turns unstable bare-D2 into a stable revival-reproducing integrator
(effect A, full revival), and (c) is the reference to validate a cheaper env against. The path to the user's
actual goal ("CTM-aware updates, dynamics more reliably, [scaling]") now runs through the **local env**:
1. **Gate the NTU {2,2} vertical patch** (NTU==exact for :up/:down as the patch→max), reusing the
   exact-cluster builder as the reference — the one missing piece before NTU is a 4-direction trajectory engine.
2. **Validate NTU against exact-cluster at D=2** on the effect-A trajectory (does NTU reproduce the
   ±2% stable revival the full-torus env gives? — a real trajectory discriminates where frozen bonds did not).
3. **Run the D=3 effect-B test with NTU** (cheap, local, scales): env-D3 vs bare-D3, both rel_floor=0.
4. (Or the true-CTM route: PEPSKit CTMRG → `bondenv_fu`, the literal "Corner Transfer Matrix" ask, χ-controlled.)
Alternatively, **bank effect-A as the deliverable** (a validated stable D=2 integrator) and scope NTU/CTM as
the next milestone. ← decision point surfaced to the user.

---

## NTU LOCAL-PATCH ENGINE (2026-06-09) — user chose NTU; vertical patch AUDITED + VALIDATED

User decision: **NTU local-patch is the production engine.** Audit + gate done.

### Audit (workflow wk0h6py1k, 6 agents) — vertical patch is CORRECT; the proposed gate was VACUOUS
Two NTU builders: `build_ntu_bondenv` (ntu_bondenv.jl, horizontal-only, validated) and
`_build_star_bond_env_ntu` (star_env.jl:256, the 4-DIR trajectory builder, hand-rolled vertical patch,
previously UNVALIDATED + missing the idcount==2 assert). Audit verdict:
- **Vertical patch is geometrically CORRECT** (3 independent lenses): both threading plaquettes present for
  up/down (worked 4×4 examples); boundary traces clean, no self-aliasing (offline idcount==2 passes all
  dirs×centers); routing/dual/assembly byte-identical to validated paths; **PSD preserved** (herm_resid=0,
  min-eig +1.2e-2 up / +3.8e-3 down on a real frozen state); **NTU-torus==exact bitwise 0.0** for up+down.
- **The proposed torus-equality gate is VACUOUS** (verify phase, both angles): on a torus patch every leg is
  internal → the boundary-trace code is DEAD → a torus gate re-validates only the routing (G1 already does
  that) and is BLIND to the P6 boundary trace. Fix demanded: a **value-level transpose oracle for vertical
  P6**. Also caught a real assert tautology in the design sketch (`>5*nsites-1`; use `>4*nsites+nsites`).

### Implemented + GATED (ntu_vertical_gate.jl) — PASS at machine precision
Added the idcount==2 safety assert (correct form) to `_build_star_bond_env_ntu` + fixed stale comments.
Built the transpose oracle (`_ntu_env_from_array` + `build_vertical_ntu_oracle`: vertical bond → transposed
lattice where it's horizontal → independent array-ncon P6 build). D=2 frozen t=1.8:
- **GV1 vertical-P6 transpose oracle: up=6.17e-14, down=1.85e-14** (≪1e-10). The vertical boundary trace is
  independently validated at VALUE level — closes the vacuous-gate gap.
- **GV2 NTU-P6 no-op-at-BIG: R=1.05e-15 U=6.11e-16 L=2.78e-16 D=3.89e-16** (all ≪1e-9). Production no-op holds.

**⇒ The cheap 4-direction NTU-P6 engine is VALIDATED.** Next: NTU-P6 D=2 trajectory (`b3p8slthf`, fork@1.6→2.8,
same fork state as the exact-cluster window) — does the LOCAL engine reproduce the exact-cluster ±2% revival
at a fraction of the cost? Then the D=3 effect-B run the full-torus env couldn't reach.

### NTU-P6 D=2 TRAJECTORY (`ntu_window_P6_dt02_to2.8.out`) — STABLE but ~14% LOW (cruder than exact-cluster)

NTU-P6 is **5.5× cheaper** (62 min vs 5.7 h, 3.5 vs 7.5 GB) and **stable** (minkept=2, nfb=0, finite, no
blowup; reproduces the revival SHAPE + timing, peak ~t=2.7) — but it **systematically undershoots**:

| t   | NTU-P6  | exact-cluster | ED      | NTU−exact | NTU−ED |
|-----|---------|---------------|---------|-----------|--------|
| 1.6 | 0.14480 | 0.14480       | 0.15558 | 0 (fork)  | −0.011 |
| 1.8 | 0.20158 | 0.22321       | 0.24370 | −0.022    | −0.042 |
| 2.0 | 0.26608 | 0.31652       | 0.33735 | −0.050    | −0.071 |
| 2.4 | 0.38531 | 0.46470       | 0.46696 | −0.079    | −0.082 |
| 2.6 | 0.41706 | 0.49275       | 0.48306 | −0.076    | **−0.066 (−14%)** |
| 2.8 | 0.41733 | 0.48253       | 0.46230 | −0.065    | −0.045 |

The NTU−exact gap GROWS to ~t=2.4 then plateaus at −0.07/−0.08: the open-boundary P6 patch is a genuinely
cruder environment whose per-step error ACCUMULATES. This is consistent with the prior frozen-bond finding
(P6-vs-exact env rel-dev 0.949, "nearly orthogonal") — the single-bond test couldn't see it, the trajectory
does. **So NTU-P6 achieves effect-A STABILITY (removes the blowup, reproduces the revival shape) but at a
~14% accuracy floor vs exact-cluster's ~2%.** On a 4×4 there is NO valid intermediate patch (P10/P12
self-alias on the wrap, need Nc≥5) to interpolate P6→exact. Candidate cheap fix: `boundary=:lambda`
(SUWeight mean-field dressing of the patch boundary — build_ntu_bondenv supports it; `_build_star_bond_env_ntu`
currently uses `:trace`). Strategic question (under review whlk8tdgd): NTU-as-cheap-stabilizer vs `:lambda`
boundary vs pivot to χ-controlled CTMRG (`bondenv_fu`, the literal CTM ask, accurate AND feasible at D≥3).

### NTU-P6 result review (2026-06-09, workflow whlk8tdgd, 3 agents) — undershoot REAL, pivot to CTMRG
- **Undershoot is REAL, not a bug** (both verify angles): identical fork/dt/alg/schedule; gap grows smoothly
  (coherent per-step rank-2 subspace-bias accumulation); nfb=0 proves the NTU env was consumed every step (no
  silent bare-SVD fallback); the SUWeight ledger lockstep is env-independent. The `:trace` open boundary closes
  with a max-mixed identity → flatter env → biases toward the bare-SVD subspace (which deterministically
  undershoots), so P6 inherits a same-signed undershoot that compounds.
- **`boundary=:lambda` REFUTED** as the fix (frozen-bond reviewer_ntu_probe data): P6-lambda rel-dev 1.234
  (WORSE than trace's 0.949), |dn| 660× worse — it double-counts the SUWeight already in the bond tensors.
- **Forward = pivot to χ-controlled CTMRG** (the literal CTM ask): χ is an env-fidelity lever independent of
  patch geometry; bondenv_fu already composes value-correct (ctm_bond_proto: fid 1.0, |dn| 1.5e-8, 26× better
  than bare, byte-identical env space to exact-cluster). `should_test_d3_ntu=false`, `should_test_lambda=false`.

### ⚠️ Subtlety flagged (not in the review): CTMRG = thermodynamic limit, NOT the finite 4×4
CTMRG (`leading_boundary` on the iPEPS unit cell) gives the **infinite-tiling** env; the ED oracle +
exact-cluster are the **finite 4×4 torus**. Different BC → CTMRG cannot be cleanly validated against the 4×4
ED the way exact-cluster was; its value is at SCALE (no oracle). The χ-ladder also measures the **finite-vs-
infinite env gap**, which bears on the open question of whether the t≈2.6 revival is a finite-size artifact.

### User decision (2026-06-09): run the CTMRG χ-ladder probe FIRST (`ctm_chi_ladder.jl`, `bv1egrsyi`)
Frozen D=2 rise bond, χ∈{8,16,24,32}: CTM bondenv_fu env vs the exact-cluster finite-torus reference.
Metrics: reldev_inf (same metric that gave NTU-P6 0.949), 1−overlap (scale-invariant), single-bond |dn|.
Read: CTM reldev_inf ≪0.949 and shrinking with χ ⇒ env is a real accuracy lever (pivot GO); a nonzero
plateau ⇒ that floor is the finite-vs-infinite BC gap (decision-critical before any multi-hour CTM trajectory).

### χ-ladder RESULTS (2026-06-09) — CTMRG VIABLE but the frozen-bond shortcut is EXHAUSTED
**Attempt 1 (t=1.8 bare, `ctm_chi_ladder_d2.out`): FAILED — CTMRG won't converge on the blown-up state.**
The bare D=2 t=1.8 state (n=0.434, pathological/near-non-injective) gave erratic reldev (0.036→0.70→1.99→1.01)
and exploding wall time (76s→1435s) — the non-convergence signature. CTMRG needs a healthy PEPS.

**Attempt 2 (t=1.0 healthy, `ctm_chi_ladder_t1.0.out`): CTMRG CONVERGES but the metric is gauge-confounded.**
CTMRG converged cleanly (err~1e-10, identical obj 2.8566e-3 across χ∈{8..48}). YET reldev_inf (1.70,1.98,1.74,
1.99,0.136) and 1-overlap (0.556,0.045,0.492,0.015,0.0092) are ERRATIC across χ — because the CTM env carries
a **gauge freedom not fixed by the objective**, so raw-tensor comparison to exact-cluster is meaningless. The
gauge-INVARIANT quantity, the truncation `|dn|_ctm`, is **rock-stable at 4.46e-8 and beats bare-SVD (3.97e-7)
by 9×** (consistent with ctm_bond_proto's 26×).

**Net:** (+) CTMRG is VIABLE — converges on healthy states, CTM-aware truncation stable & better than bare.
(−) The frozen-bond probe CANNOT settle CTM accuracy: env-comparison is gauge-confounded; single-bond |dn| is
non-discriminating (near-lossless); the healthy-high-entanglement regime is D≥3-only (exact-cluster infeasible
there); and CTMRG models the infinite-tiling, not the finite 4×4 ED. **The cheap probes are exhausted — only a
CTM-aware TRAJECTORY vs ED can settle accuracy**, and that is a real build (amortize the CTMRG env rebuild
once per sweep; CTMRG was 75s–1461s/build at D=2 here, so per-bond-per-star is infeasible, per-sweep is ~ok).

### Strategic synthesis (surfaced to user 2026-06-09)
The 4×4 testbed has DELIVERED its science: a loop env makes D=2 stable+~2% (effect A); the lever is D (dt-
converged); the env-engine cost/accuracy landscape is mapped (exact-cluster accurate-but-D=2-only; NTU-P6
cheap-but-~14%; CTMRG viable, scalable, thermodynamic-limit). Further 4×4 frozen-bond work hits diminishing
returns. Fork: (A) bank the 4×4 story + scope CTMRG-at-scale as future work; (B) build the amortized CTM-aware
trajectory engine (thermodynamic limit; validated by χ-convergence + cross-checks, not the finite 4×4 ED).
