# CTM-Upgrade Execution Workflow — review → improve → validate (full-backlog edition)

> Status: workflow design (2026-06-07). This is the **orchestration layer** — it sequences
> and gates the work; it does **not** restate the physics. Read alongside the technical docs:
> - `docs/superpowers/plans/2026-06-04-ctm-aware-evolution.md` — the PR-0→PR-5 / P0–P6 / G0–G7 plan.
> - `docs/superpowers/specs/2026-06-04-pepskit-native-ctm-aware-migration-design.md` — Stage A/B/C/D.
> - `notes/stage2-truncation/2026-06-04-ctm-aware-evolution-state.md` — engine fork + prototype state.
> - `notes/2026-06-02-package-review-and-roadmap.md` — the backlog this clears.
> - `notes/methodology/revival-validation.md` — the validation law (load-bearing).
>
> **Decisions locked for this workflow (2026-06-07, author):**
> 1. **Scope = full review backlog** — CTM spine **+** heaviness slim **+** test de-trapping **+** Stage-3.
> 2. **First environment engine = exact finite-cluster** (state-doc §6 recommendation), with
>    **NTU finite-patch as the production follow-on**. Per-star CTMRG stays a diagnostic ceiling only.

> ## ⚠️ SUPERSEDED PENDING DE-RISKING SPIKE (2026-06-07 fan-out review)
> A 69-agent fan-out review (10 verified blockers) found load-bearing errors below. See
> **`notes/stage2-truncation/2026-06-07-workflow-fanout-review-verdict.md`** for the full verdict.
> Engine decision #1 is **UN-LOCKED** pending an N4+N2+controls spike. Corrections that hold regardless:
> - "reuses the validated oracle contraction / zero new infrastructure / runnable now" (§4 PR-0/PR-3,
>   §12) is **false** — the exact-cluster `{2,2}` env is a heavier **double-layer** object; the
>   single-layer oracle was killed at D=4 (`g0_frozen_step.out`).
> - **G1** (§4/§5) is wired into `test_d_convergence.jl` = the **legacy** bare-SVD backend, not the
>   native path PR-3 repairs. Re-author onto the native backend; first prove the native oracle survives
>   a full D=4 trajectory.
> - **A1's monotone D-ladder is legacy-only**; the native ladder is **flat at D≥3** and no dt-control
>   exists. Add a dt-convergence gate before crediting "environmental, not Trotter."
> - **`exact_pxp_energy_density_finite` already exists** (`Observables.jl:926`, exported) — G3/A8 needs a
>   `::PXPIPEPSState` method, not green-field.
> - A2's "needs no `positive_approx`" is false (`fixgauge_benv` consumes `positive_approx`'s `Z`);
>   keep it + a cond(Z'Z) tripwire. Decouple PR-1 (CTM η is "CTM-FU path only") from PR-3.
> - The design-doc names **NTU the primary engine** (`...evolution.md:41`, PR-3=NTU); this doc inverted it.

---

## 0. What "done" means (exit criteria)

From the migration spec's exit criteria + the review's three stages:

- [ ] **CTM-aware milestone** — one env-aware truncation shows a measured revival-rise-lag
  reduction at fixed D (G1) **and** fixes the native D=2 blowup, judged by the **exact-oracle
  trajectory RMS over [0, 2.8]**.
- [ ] **Native D>1 validated** vs the exact oracle by the trajectory metric (largely landed in
  `b33b416`; close the residual Stage-B items).
- [ ] **Trust is real** — CTM gates on a true η residual, not the `truncation_error` mislabel;
  no test pins the D>1 simple breakdown as "expected."
- [ ] **Slimmer** — TFIM axis + duplicate drivers + dead exports removed/relocated per the
  delete-vs-demote call; suite still green.
- [ ] **Stage-3 real** — a return-amplitude/Loschmidt revival observable + a search layer ranked
  vs the Néel baseline; the `:z_up`/seam crash is legible.
- [ ] **Single engine** — ITensors/Strided dropped, dual-engine selector removed (final, gated on
  D>1 parity).

---

## 1. The meta-loop — "continuing the review process," made repeatable

Every slice below runs the same four beats. This *is* the review process, instantiated per PR:

1. **REVIEW** — scope the slice. For *research-grade* claims (cost model, env geometry,
   metric well-posedness) fan out independent agents and **adversarially verify before
   committing** (the 19-agent pattern that produced the plan). For *code* slices, a focused
   diff read (`/code-review` style).
2. **IMPROVE** — TDD: write the **red physics gate first**, then the **minimal** impl, always
   **behind an opt-in flag**. Never touch the legacy ITensors backend or the bare-SVD default.
3. **VALIDATE** — fast gate every pass; the milestone's physics gate (G0–G7); nightly for the
   full suite.
4. **RECORD** — audit note to `notes/<stage>/`; **surface owner forks, don't decide them**.

### The validation law (non-negotiable — `notes/methodology/revival-validation.md`)

> - Headline metric = **`exact_density_finite` trajectory RMS over [0, 2.8]**. Never CTM
>   (χ=8 *flatters* the revival by ~the signal size). Never a single time (t=2.6 is a
>   curve-crossing that inverts the ladder). **D=1 is never validation.**
> - Compare **observables, never raw tensors**.
> - **Self-deception guards**: (a) the instrument that drives truncation (the env) and the one
>   that judges (the exact oracle) must differ; (b) energy-bounded **and** density-improved
>   *jointly* (a too-aggressive clamp looks energy-stable while damping the revival);
>   (c) energy-conservation / reverse-echo are necessary-not-sufficient — the exact-oracle
>   comparison is the only sufficiency check.

---

## 2. Track 0 — the iterate harness (build FIRST, ~1 slice)

The review found the suite was *unloadable* until 2026-06-02; the gate must be fast enough to
run every pass and trustworthy enough to authorize a merge.

- **Fast gate (per pass, target seconds in a warm process):**
  `test_pxp_validation.jl` + `test_pxp_d2_localization.jl` + `test_star_simple_update.jl`.
  Run via the existing ARGS subset:
  `julia --project=. test/runtests.jl test_pxp_validation.jl test_pxp_d2_localization.jl test_star_simple_update.jl`.
  Kill the ~75 s compile tax with a **warm Julia daemon / sysimage** (Revise-style).
  → *Refinement of the review's gate:* `test_ipeps_compression.jl` (582 s) is too slow for a
  per-pass gate; demote it to the **pre-merge** gate, not the iterate gate.
- **Pre-merge gate:** the fast gate + `test_ipeps_compression.jl` + `test_d_convergence.jl`
  (the trajectory-RMS D-ladder) + whichever physics gate the slice touches.
- **Nightly:** full suite + `SQUAREPXP_EXTENDED_TESTS=1` + Aqua + the **real-CTM-vs-ED D>1**
  test (Track A) + long-time D=4 conditioning. Under load (`uptime ≳ 100`) use
  `JULIA_TEST_NWORKERS=8` and capture `PIPESTATUS` (the pmap harness under-reports failures).

**Exit:** fast gate runs warm in < ~30 s; nightly is wired in CI (extend `.github/workflows`).

---

## 3. Track map & dependency graph

```
            ┌─────────────────────────────────────────────────────────────┐
 Track 0 →  │  iterate harness (fast / pre-merge / nightly)                │  enables all
            └─────────────────────────────────────────────────────────────┘
                 │
   ┌─────────────┴───────────────┐
   ▼                             ▼
 Track A (trust + de-trap)     ── concurrent prereq ──┐
   PR-1 real η                                        │
   de-trap golden tests                               │
   real-CTM-vs-ED D>1 test                            │
   promote trajectory-RMS to prod aggregation         │
   ┌───────────────────────────────────────────────── │ ───────────────┐
   ▼                                                   ▼                 │
 ★ CRITICAL PATH — Track B (CTM-aware evolution)                        │
   PR-0/G0  decision experiment  ──┐ (needs nothing but the oracle)     │
   P0       doc PSD fix (blocks env merges)                             │
   PR-2     shims + convention-B write-back + P6 tripwire ◄── Track A   │
   PR-3 ★   EXACT FINITE-CLUSTER env → star split (opt-in) ◄── PR-0,P0,PR-2
   PR-4     NTU finite-patch (production) ◄── PR-3 as gold standard     │
   PR-5     Stage-D 5-tensor ALS (gated on P4/P5 + a measured plateau)  │
   └────────────────────────────────────────────────────────────────────┘
   ▲                                                   ▲
   │ shared building blocks (build once, §6):          │
   │   exact_pxp_energy_density_finite (G3)            │
   │   Loschmidt / reverse-echo observable (G4 + Stage-3 objective)
   │
 Track C (heaviness slim)  ── parallel fan-out, gated only by "suite stays green"
 Track D (Stage-3 ScarFinder) ── parallel; reuses the Loschmidt observable from §6
```

Critical path = **Track 0 → PR-0 → P0 → PR-2 → PR-3**. Everything else parallelizes against it.

---

## 4. ★ Critical path — Track B (CTM-aware evolution)

Each PR: **red gate → minimal impl behind an opt-in flag → audit note**. Legacy retained throughout.

| PR | Red gate (write first) | Minimal impl — hook | Exit gate | Owner checkpoint |
|---|---|---|---|---|
| **PR-0** `G0` decision experiment | assert (i) exact-untruncated, (ii) exact-best-rank-D, (iii) simple-update densities, report (ii)−(i), (iii)−(i) at t∈{1.6,1.8,2.0,2.2}, D=2,4 | `scripts/g0/` harness over `exact_density_finite`/`_native_dense_state` (`PEPSKitObservables.jl:474-553`) + the exact star gate. **Also persist** the ephemeral `/tmp/probe{2,4}.jl` 95–350× ratios into a committed artifact. | branch decided (see note ▼) | — |
| **P0** doc PSD fix | n/a (doc) | Correct the spec's "indefinite N ill-posed" (design doc :177-187) **or** roadmap :66-67; one-line: the truncation objective is the **PSD norm-fidelity** (`cost_function_als`/`bondenv_fu`). | **no env code merges until one doc is fixed** | — |
| **PR-2** trust shims + writeback | per-primitive contract test (symmetric √S split) for each PEPSKit private used; convention-B equivalence test under a non-Schmidt write-back; **P6** D>1 trajectory blockade tripwire | one internal shim module for `absorb_s`/`flip_svd`/`bondenv_fu`/`bond_truncate`/`_qr_bond`/`positive_approx`/`calc_convergence`, pin `PEPSKit=0.7.0`; local-SVD diagnostic-ledger re-extraction; relabel `measure_simple` "tree diagnostic" and **exclude from every gate** | shims fail loudly on a PEPSKit bump; convention-B equivalence holds to 1e-10 | P3: keep `SUWeight` as diagnostic-only ledger vs deprecate — **default keep** |
| **PR-3 ★** exact finite-cluster env | **no-op at maxdim=BIG** reproduces the lossless density; env vs brute-force `ncon` at D=1,2; then **full D=2 & D=4 trajectory** vs ED/legacy | Build the `{2,2}` `BondEnv` by the **same finite-torus double-layer contraction as `_native_dense_state`**, two bond-end virtual legs left open, physical legs traced → **PSD by construction**. Feed `PEPSKit.bond_truncate(...; ALSTruncation)` via the §5 state-doc idiom + `fixgauge_benv` (cond(Z'Z)~3.4e20 on the rise). Replace the bare `_split_bond` at **`PEPSKitStarUpdate.jl:322`** behind an **opt-in flag** on `project_star_pepskit!`/`evolve_pepskit!` (driver `PEPSKitEvolution.jl:143-189` has no env hook today). Rebuild `N_dir` per arm (peel staleness, §4a). | **G1** (fixed-D lag ↓, D-monotone) + **G2** (monotone in patch radius, plateaus) + **G5** (holds on 4×4 *and* 6×6) + **G7** (sublattice Z₂) — **fixes the D=2 blowup as a bonus** | — |
| **PR-4** NTU production engine | same G1/G2/G5; **cost-gate vs PR-3** (accuracy-per-cost) | hand-contract a finite plaquette-ring `{2,2}` patch env (PSD, no CTMRG, no `rotl90` pain); same hook/flag as PR-3 | passes on **6×6** where exact-cluster is infeasible | engine for cells > 4×4: NTU vs warm-CTMRG — **default NTU** |
| **PR-5** Stage-D ALS | **G6** snake-MPO path-independence; joint-vs-per-bond optimum differs measurably; **G1 below the PR-3/PR-4 plateau** | wire M7 (`PEPSKitStarSnake.jl:78-106`) to *apply* the gate; author `E₊` (no PEPSKit primitive); 5-core ALS | beats PR-3/PR-4 at the entangled peak | **gated**: only if M3/M4 leave a plateau worth P4's 10–50× cost |

▼ **PR-0 honesty note (reconciling the plan's G0 with state-doc §4).** G0 settles **staging**
(is there headroom in a better *truncation metric*: does exact-best-rank-D (ii) beat simple-update
(iii)?). It does **not** demonstrate the *fix*: the D=2 blowup is **cumulative** (a single frozen
step's truncation error is ~5.4e-4 even when the frozen state is already corrupted to n=0.434 vs
ED 0.156, `g0_frozen_step.out`). The fix is demonstrated **only** by PR-3's full-trajectory run
from t=0. So: run G0 to choose the branch and to commit the load-bearing probe numbers; treat
PR-3's trajectory as the actual proof. §2's evidence predicts the "(ii)≈(i), (iii) lags both"
branch → the error is the **environment**, not the bond budget → exact-cluster/NTU/Stage-D, not a
bare reweighted SVD.

**Why the chosen engine de-risks PR-3:** the exact finite-cluster env is built by a *direct finite
contraction*, not a CTMRG fixed point — so it **sidesteps the vertical-bond `rotl90(env)` blocker**
that stalled the per-star CTMRG prototype (`ctm_star_truncate.jl`), needs no `positive_approx`
(PSD by construction), and reuses the already-validated oracle contraction. It is the **gold
standard on 4×4** to prove the thesis; NTU (PR-4) is the cheaper approximation that scales to 6×6.

---

## 5. Track A — trust plumbing & de-trapping (concurrent prereq, gates B's credibility)

Runs alongside PR-0; **must land before PR-3's validation is trustworthy**.

1. **PR-1 — real CTM η (P1).** `leading_boundary` info has no η/converged field; the adapter
   aliases `truncation_error`→`residual` (`PEPSKitMeasurements.jl:493,497-499`) and the native
   context never inspects convergence (`PEPSKitObservables.jl:181,199-210`). Compute
   `η = calc_convergence(...)` (`ctmrg.jl:191-207`); store η / truncation_error / condition_number
   as **three distinct fields**; gate `converged := η ≤ tol`. **Rewrite the tests that assert the
   mislabel** (`test_pepskit_measurements.jl:167,340`; `test_ctm_trust.jl:286-318`) or the bug is
   frozen by the suite. *Red gate:* η matches the loop on a converged env, large on `maxiter=2`.
2. **De-trap the golden tests.** The `>1e-4` lower-bound + `@test_broken` tests **pin the D>1
   simple breakdown as expected** and trap any fix. Convert to `@test_broken` + add
   `exact_finite < 1e-6` **upper-bound** gates.
3. **Add the missing real-CTM-vs-ED D>1 test** (every CTM-value test today uses fake closures).
4. **Promote trajectory-RMS into production aggregation.** `validate_pxp_convergence` aggregates a
   single-time **max**; the trajectory metric lives only in `test_d_convergence.jl:67-122`. Lift it
   into the production validator (G7 / G5 reuse it).

---

## 6. Shared building blocks (build once — used by Track B **and** Track D)

These appear as "BUILD" items in multiple gates; author them once, early (in/after PR-2):

- **`exact_pxp_energy_density_finite(::PXPIPEPSState)`** — mirror `_native_dense_state`, apply the
  star H per center, ≤16 sites, ledger-independent. Powers **G3** (energy-drift tripwire) and the
  Stage-3 energy-variance proxy. *Does not exist today* (native energy paths are mean-field or
  CTM-contaminated).
- **Loschmidt / reverse-echo / return-amplitude observable** — a gauge-invariant exact-contraction
  distance (4×4) or frozen-ED endpoint (6×6). Powers **G4** (ED-free reverse-echo for the CTM
  track) **and** Stage-3's real revival objective (today `RevivalObjective` is just instantaneous
  `|n_even − n_odd|`, `ScarFinder.jl:544-552`). One observable, two consumers.
- **PEPSKit-private shim module** (from PR-2) — the single import surface for env work.

---

## 7. Track C — heaviness slim (~3–4k lines; parallel fan-out, gated by "suite green")

Low-risk; ideal for **agent fan-out** (one agent per target: delete/move → fast gate → report).
Sequence so it never blocks the critical path. **Coordinate with Track D**: the ScarFinder-objective
deletion is *deferred* — those placeholders get *implemented* in Stage-3, not deleted.

| Target | ~Lines | Action | Wave |
|---|---|---|---|
| TFIM axis (`FiniteTFIM*`, TFIM model/obs/benchmarks) | ~950 src +700 test | move to `examples/` or delete (orthogonal to PXP) | 1 |
| Dead exports + 245-symbol surface | ~30–40 exports | strip; relax `test_public_docs` to core API (loop-until-dry census) | 1 |
| 30 design docs (kagome/PESS never built) | — | archive to `notes/archive/` | 1 |
| Duplicate PXP campaign drivers + struct families + serialization | ~550 | merge into one grid runner | 2 |
| `*_simple` individual exports (D>1-divergent) | ~15 exports | demote from public API | 2 |
| `LowVarianceObjective`/`TargetEnergy` placeholders (`ScarFinder.jl:285`) | ~150 | **defer** → implemented by Track D, not deleted | — |

**Owner checkpoint:** aggressive **delete** vs conservative **demote-to-experimental**. *Default:*
delete TFIM + dead exports (wave 1), demote `*_simple` (wave 2). The "662-line `fix_bond_gauge!`
orphan" in the review is **already deleted** — skip it.

---

## 8. Track D — Stage-3 ScarFinder (parallel; depends on §6 Loschmidt observable)

The review's verdict: the search is *not expressible today* — `scarfinder!` is a single-seed
rank loop with no optimizer, no Néel baseline, and no real revival metric. Path to real Stage-3:

1. **Legible crash** for `:z_up` / odd-cell checkerboard seams (blockade-forbidden ‖θ‖≈0 currently
   throws opaquely deep in the star split) — clear error or graceful skip.
2. **Real revival objective** — the §6 Loschmidt/return-amplitude observable over the trajectory
   time series (replaces instantaneous staggered-mag, which Néel maximizes at t=0).
3. **`scarfinder_search` outer layer** over a candidate-state family, ranked **vs the Néel
   baseline**.
4. **(optional) energy-variance proxy** (the standard scar selector) via §6's exact energy density.
5. **Bug fixes:** `ScarFinderAudit.jl:338` (`product_square_ipeps` → `_build_state`, silently drops
   reversibility for all checkerboard rows); `scripts/run_scarfinder_audit.jl` default
   `initial_state=:z_up` (driver broken out of the box).

**Owner checkpoints (forced when Track D starts):**
- **"Better than Néel" metric** — longer-lived revival vs lower energy variance vs larger
  return amplitude. *Default:* return-amplitude (reuses §6, the cleanest "revival").
- **Search space** — per-sublattice product-angle family vs low-D seeds. *Default:* product-angle.
- **Geometry** — even-Lx/Ly + `:serial` (Néel works today) vs `:five_color` (can't host Néel).
  *Default:* even-cell `:serial`.

---

## 9. Owner decision checkpoints (one table; defer each until its PR is reached)

| # | Decision | Forced at | Default |
|---|---|---|---|
| 1 | Engine #1 | — | **exact finite-cluster** ✅ (locked 2026-06-07) |
| 2 | Scope | — | **full backlog** ✅ (locked 2026-06-07) |
| 3 | PSD-doc target (spec vs roadmap) | P0 | fix the **spec** :177-187 |
| 4 | `SUWeight` ledger: keep-as-diagnostic vs deprecate | PR-2 (P3) | **keep** as diagnostic-only |
| 5 | P·U vs P·U·P under env truncation (P6) | PR-3 | start **P·U**, tripwire-gate; escalate to P·U·P only if blockade leaks |
| 6 | Stage-D go/no-go | after PR-3/PR-4 | **no-go unless** a measured plateau > P4 cost |
| 7 | Slim: delete vs demote | Track C | delete wave-1, demote wave-2 |
| 8 | "Better than Néel" metric / search space / geometry | Track D | return-amplitude / product-angle / even-`:serial` |
| 9 | Production engine for cells > 4×4 | PR-4 | **NTU** |

---

## 10. Multi-agent orchestration map (where fan-out pays)

- **Solo, focused:** PR-0/G0 (one experiment), PR-3 hook (one surgical replacement), P0 (doc).
- **Adversarial verification (19-agent style — N independent skeptics, majority-refute kills):**
  the research-grade claims — **P4** cost model (the brainstorm's "single biggest missing input"),
  **P5** `E₊` geometry/aliasing, the indefinite-metric spectrum probe, and the PR-0 branch
  interpretation. These decide expensive staging; verify before committing.
- **Parallel fan-out (pipeline: transform → fast-gate → report):** Track C slim (one agent per
  deletion target, loop-until-dry on the dead-export census) and Track A de-trapping (one agent per
  test file).
- **Contract-test pass:** PR-2 shims and the PR-3 env-construction correctness (no-op-at-BIG +
  vs-`ncon` at D=1,2) — a verifier per primitive.

The harness for all of the above is the repo's existing workflow infrastructure
(`scripts/run_autonomous_loop.sh`, the superpowers pattern) plus the fast/pre-merge/nightly gates.

---

## 11. Milestone sequence (headline = M3)

| Milestone | Contains | "Done" signal |
|---|---|---|
| **M0** enabling | Track 0 harness; **P0** doc fix; persist `/tmp` probe artifacts | fast gate warm < 30 s; one PSD doc corrected |
| **M1** decision | **PR-0/G0**; ∥ **PR-1** real η + de-trap | C-vs-D branch chosen; no test pins D>1 breakdown |
| **M2** trust+blocks | **PR-2** shims + convention-B writeback + P6 tripwire; §6 building blocks; ∥ Track C wave-1 | convention-B equivalence 1e-10; suite green minus TFIM |
| **M3 ★** the upgrade | **PR-3** exact-cluster env → star split (opt-in); full D2/D4 trajectory validation | **G1+G2+G5+G7 pass; D=2 blowup fixed; D=4 lag ↓** |
| **M4** production | **PR-4** NTU, cost-gated, on 6×6; ∥ Track D revival observable + search + legible crash | NTU matches exact-cluster within cost budget on 6×6 |
| **M5** (gated) | **PR-5** Stage-D ALS | only if M3/M4 leave a plateau worth P4 cost |
| **M6** final | drop ITensors/Strided, remove dual-engine selector; Track C wave-2 | single engine; D>1 parity held |

---

## 12. Start here (next concrete action)

**PR-0 / G0 is runnable now** — zero new infrastructure, reuses the validated oracle, and it is the
gate that decides the entire C-vs-D staging. Concretely:

1. Promote `scripts/g0/frozen_step_diagnostic.jl` into a committed G0 harness that, at
   t∈{1.6,1.8,2.0,2.2} and D∈{2,4} on the 4×4 Néel cell, reports densities (i) exact-untruncated
   post-gate, (ii) exact-best-rank-D, (iii) simple-update — and emits `(ii)−(i)` and `(iii)−(i)`.
2. Persist the `/tmp/probe{2,4}.jl` 95–350× error/discard ratios into
   `notes/stage2-truncation/data/` (currently load-bearing evidence that lives only in `/tmp`).
3. In parallel (independent), open **PR-1** (real η) and the de-trap pass — they unblock PR-3's
   credibility and touch disjoint files.

Author sign-off needed only on decisions #3–#9 as each PR reaches them; #1 and #2 are locked.
