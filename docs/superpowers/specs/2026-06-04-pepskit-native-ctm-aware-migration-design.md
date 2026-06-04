# PEPSKit-native + CTM-aware PXP dynamics: migration design

Date: 2026-06-04. Author: architecture review (13-agent SuperPower scan + adversarial
physics diagnosis + first-hand source/PEPSKit-internals verification).

Extends — does **not** replace — `docs/pepskit-native-refactor.md` (the existing M0–M7
plan) and the Stage-2 backlog in `notes/stage2-truncation/`. This doc is the skeptical
audit of that plan plus the staged path to CTM-aware evolution.

## Goal

Preserve the genuine five-site square-PXP **star** physics while (1) making PEPSKit
`InfinitePEPS` the single source-of-truth state, (2) demoting `SUWeight`/λ to an update
aid + mean-field diagnostic, (3) making CTMRG the source of physical norms/observables,
and (4) reaching **CTM-aware / full-environment evolution** (the real fix for the
revival-rise lag). Do this incrementally, behind physics validation gates, with the
legacy ITensors backend retained until D>1 parity is proven.

## Repository scan summary (what is what, verified)

**PEPSKit-backed (native, branch `pepskit-native-backend`):**
- `PEPSKitBackend.jl` — `PXPIPEPSState = InfinitePEPS + SUWeight` + bookkeeping
  (`mutation_version`, `log_norm`); product/checkerboard constructors; coordinate maps.
  Convention "B": every stored Γ carries √λ, so `measurement_peps(state) === state.peps`
  (`:301`). **InfinitePEPS is already the source of truth.**
- `PEPSKitStarUpdate.jl` — native five-site star gate (`pxp_star_gate_tensormap`) and
  `project_star_pepskit!` simple update. Uses PEPSKit `absorb_weight` + TensorKit
  `tsvd/leftorth/truncdim/truncbelow`, and **unexported** `PEPSKit.absorb_s`,
  `PEPSKit.flip_svd` (`:228,:230`).
- `PEPSKitEvolution.jl` — `evolve_pepskit!` Trotter driver, five-colour palindromic
  schedule, deliberately **not** PEPSKit generic `time_evolve`; `reverse` flag (real-time).
- `PEPSKitObservables.jl` — two layers: (i) CTM observables via `LocalOperator` +
  `leading_boundary` + `expectation_value` (norms from CTMRG); (ii) CTMRG-free
  simple/local observables on convention-B patches + the `measure_simple`/`evolve!`/
  `reverse_evolve!` adapters + `pxpipeps_from_square_ipeps` converter.
- `PEPSKitStarSnake.jl` — M7 exact (rank≤2) MPO factorization of the star gate; **not
  wired into evolution** (correctly fenced off).

**Still custom ITensors / Γ–λ (the production path, default `:legacy_itensors`):**
- `SquareIPEPS.jl` (Γ–λ state, **full-λ** absorb/deabsorb), `StarSimpleUpdate.jl`
  (`project_star!`, the full-λ deabsorb hazard, `canonicalize_simple!`),
  `IPEPSEvolution.jl` (legacy Trotter), `StarModels` in `PXPModel.jl` (ITensors gate
  wrappers), the simple-observable half of `Observables.jl`.
- `PEPSKitMeasurements.jl` (1267 lines): ~150 conversion glue (`to_pepskit_infinitepeps`),
  ~155 operator→`LocalOperator` builders (stay needed), ~860 genuine CTM logic.

**Dead / transitional (census, corrected against the refactor doc):**
- `CTMGaugeReadiness.jl` / `fix_bond_gauge!` are **already fully deleted** (grep/`find`
  return nothing; only a stale comment `SquareIPEPS.jl:281` + stale notes remain). The
  "662-line orphan to wire-or-delete" item in `notes/2026-06-02-package-review` is **stale**.
- `canonicalize_simple!` — orphaned from `evolve!` **on purpose** (kept as a tested Vidal
  primitive for diagnostics + the future full-update); empirically ruled out as a
  dynamics fix (`notes/.../2026-06-02-stage2-regauge-map.md`). Do **not** re-propose it.
- `LowVarianceObjective` — genuine dead placeholder (`ScarFinder.jl:285`, never scored).
- `RevivalObjective` — **not a revival**: it is instantaneous staggered magnetization
  `|n_even − n_odd|` at one time (`ScarFinder.jl:544-552`); no Loschmidt/return amplitude
  exists in ScarFinder. (Stage-3 concern, out of scope here.)
- `to_pepskit_infinitepeps` is **NOT** "dead once states are native" (refactor-doc claim
  overstated): `measure_ctm` — the only production CTM path, used by `TrustedCTMBackend`
  and all `PXPValidation` CTM — still routes through it (`PEPSKitMeasurements.jl:808`).
  It dies only when a native CTM path replaces `measure_ctm`.

## Physics / algorithm diagnosis

**Is the simple-update dynamics trustworthy?** Yes, as a *correct* Bethe/mean-field
simple update — not a buggy one — within a clearly bounded regime:
- Collapse (t<1.4), D=2,3,4: ED-validated, near-exact.
- Revival rise (t≈1.4–2.6): a controlled **mean-field-environment lag** that shrinks
  monotonically with D (trajectory RMS D2 2.6e-2 → D4 9.8e-3). This is the dominant
  physical error and is **structural** (the single-site λ² environment is a tree/Bethe
  approximation that discards the loop correlations rebuilding Néel order).
- The t=2.6 peak is a curve-crossing; judge by the **n(t) trajectory RMS over [0,2.8]**,
  never one time. D=1 is **never** validation.

**Main failure modes:**
1. **Mean-field λ² environment** (the revival-rise lag). No λ/rel_floor/regauge knob can
   fix it — only an environment-aware truncation can. *This is why CTM-aware evolution is
   the real target.*
2. **Full-λ deabsorb poisoning** — **legacy only.** `deabsorb_link_weight`
   (`SquareIPEPS.jl:402-414`) inverts the full λ with a fixed 1e-14 floor; retained
   λ≈1e-6..1e-12 invert to ≈1e6..1e12 on the next neighbouring star. `rel_floor=1e-3`
   *conditions* (caps cond at 1/rel_floor) but does not eliminate it; the `min_lambda`
   diagnostics that measure it never gate the update.
   **Correction (verified in PEPSKit source):** the **native** backend does **not**
   inherit this. `absorb_weight` uses `pow = inv ? -1/2 : 1/2`
   (`PEPSKit/.../suweight.jl:240`) — it absorbs/de-absorbs **√λ**, with a *relative* floor
   `eps^(3/4)·λ_max ≈ 1.8e-12·λ_max` (`util.jl:48,55`). So the native update only ever
   inverts √λ (≈1e-6, not 1e12) and zeros sub-threshold weights. **The native backend
   already implements the principled Vidal √λ deabsorb the legacy notes call
   "unimplemented."** This is a physics argument *for* the migration.
3. **CTM convergence is not actually measured.** In PEPSKit 0.7.0 `leading_boundary`'s
   `info` is the last per-iteration projector info (`truncation_error`, `condition_number`)
   — **no** η residual / iteration count / `converged` flag. The adapter latches
   `truncation_error` as "residual" (`PEPSKitMeasurements.jl:490-499`), so `accepted` (and
   `tight_ctm_trust_policy`'s 1e-6 residual gate) reduce to an **SVD-discard** test, not a
   fixed-point-convergence test. Trust gates almost entirely on finite-χ **drift**.
4. **χ=8 flatters the revival** by 2.5e-3..1.3e-2 — the *same order as the signal* — so the
   headline benchmark routes around CTM to the exact 16-site oracle. A CTM-aware truncation
   that uses χ=8 would optimize against a contaminated environment.
5. **Blockade-forbidden patches** (`:z_up`, odd-cell checkerboard seams) annihilate a star
   (‖θ‖≈0); the SU aborts with a named-center error rather than continuing.

**Why CTM-aware evolution is the real target:** the revival-rise lag *is* accumulated
wrong-truncation by the mean-field environment. Reweighting the star SVD by the real
CTMRG bond environment N (so kept Schmidt vectors minimize ‖θ − θ_trunc‖_N against the
actual network) is the only thing that attacks it. The empirically-decisive experiment
(regauging is a pure gauge transform, cannot change what the next step truncates) already
ruled out the cheaper alternative.

**The five-site cross is genuine and must stay custom** (verified bit-exact:
`kron(X, P↓,P↓,P↓,P↓)`, 2 nonzeros, site order center,right,up,left,down; basis 1 =
excited; `[P,U]=0`; coordinate round-trip and five-colour disjointness numerically
confirmed, min L1 = 3). PEPSKit's generic `time_evolve` supports only 1/2/L-3-site
clusters; it cannot represent the cross without a path-biased snake. A Stage-D full update
is therefore a **5-tensor patch ALS against the whole plus-star environment**, not a
2-site bond environment.

## Target architecture (decisions)

1. **Source of truth = PEPSKit `InfinitePEPS`.** *Yes* — and it is a conditioning win
   (√λ absorb/deabsorb, relatively floored), not just maintenance. Conditional on closing
   the D>1 validation gap below.
2. **λ / `SUWeight` = strictly auxiliary** (update aid + explicitly-labelled mean-field
   diagnostic). *Not yet true in code:* the native simple-observable layer reads
   `state.weights.data` and reweights cut bonds to λ² (`PEPSKitObservables.jl:304-333`),
   and `weight_entropy` presents a Bethe quantity as a Schmidt entropy. Making λ truly
   auxiliary requires the native CTM path to become the trusted observable.
3. **CTMRG = source of physical norms/observables** — realized only in the CTM layer,
   which is currently gated *out* of the native production path.
4. **Stays custom:** the five-site star gate + `project_star_pepskit!`; the five-colour
   schedule; PXP observables; ED/exact-finite oracle. **Delegated to PEPSKit:**
   `absorb_weight` (√λ), `InfinitePEPS`/`SUWeight`, `CTMRGEnv`/`leading_boundary`/
   `expectation_value`/`LocalOperator`, `correlator`/`correlation_length`.
5. **The three endgame goals collapse into ONE work item:** *build & wire the native CTM
   measurement path for `PXPIPEPSState`.* It simultaneously (a) makes λ auxiliary,
   (b) unblocks dropping `to_pepskit_infinitepeps`/ITensors/Strided, and (c) provides the
   gauge-invariant oracle that closes the D>1 validation gap.

## Staged migration

The user's A/B/C/D maps onto reality as follows (much of A/B already exists on this branch).

### Stage A — cleanup / honesty (low risk, no new physics)
**Deliver:** (a) write the architecture decision into `notes/methodology/` incl. the
verified √λ-absorb fact; (b) add the CTM-layer Hermiticity guard (`_real_expectation`
pattern) and a **bond-dimension guard** (clear error, not silent D^12 hang) to the simple
layer `_star_energy/_local_density/_nn_excited`; (c) shim `PEPSKit.absorb_s/flip_svd/
truncrank` behind one internal module with a load-time contract `@assert` (symmetric √S
split); (d) fix the `|1><1|` density comment typo (`PEPSKitObservables.jl:58`) and rename
`weight_entropy`→ document as mean-field; (e) add a `min_lambda` abort gate to **legacy**
`project_star!` (and native), wiring the existing `touched_min_lambda` diagnostic.
**Validation:** all existing gates green; contract test fails loudly on PEPSKit upgrade.
**Risk:** trivial.

### Stage B — PEPSKit-native simple star update as the *trusted, validated* path
*(the update itself is done; this stage makes it trustworthy)*
**Deliver:** (a) a native CTM measurement path for `PXPIPEPSState` (wire the existing
`pepskit_native_context`/`measure_ctm_pepskit` into `TrustedCTMBackend` and
`PXPValidation` for native states); (b) give `exact_density_finite` (and the
`exact_*_finite` family) a `PXPIPEPSState` method by converting to per-site dense arrays
(the simple layer already reconstructs these) + reuse the `:boundary` contractor;
(c) **close the D>1 validation gap** — native-vs-ED trajectory RMS at D=2,3,4 via the
exact oracle, a convention-B equivalence test, an anisotropic-D=2 single-star
round-trip-vs-dense test, and a coordinate-reflection trap (see Validation Plan).
**Gate before any public-API migration:** native D=2 short-time `density_error_exact_finite
< 1e-3` at t=0.1 *and* the 4×4 Néel revival trajectory RMS matching the legacy ladder.
**Risk:** medium — touches the converter and the oracle typing; the √λ-vs-full-λ
convention difference is a silent-physics trap (compare observables, never raw tensors).

### Stage C — CTM-aware bond truncation (first meaningful CTM-aware evolution milestone)
**Deliver:** reweight the native star-bond split (`PEPSKitStarUpdate.jl:_split_bond`,
`:207-232`; legacy hooks `StarSimpleUpdate.jl:334/342`) by a bond environment N built from
a **converged, large-χ** `CTMRGEnv`, replacing the λ² mean-field metric.
**Hard prerequisites (do NOT code before these):**
- Fix CTM residual capture (Stage A/validation) — N must come from a *verified-converged*
  environment, χ large enough that environment error ≪ truncation error.
- **Resolve the indefinite real-time metric objective** — in real time N is
  non-Hermitian/indefinite, so ‖θ − θ_trunc‖_N is ill-posed; choose explicitly between
  regularized eigh + pseudo-inverse on near-null directions vs an imaginary-time/PSD
  surrogate (owner decision, `notes/.../pxp-improvement-brainstorm.md`).
- **Convention-B reconciliation** — post-CTM-truncation bond weights are *not* a Schmidt
  spectrum; define whether the `SUWeight` ledger is re-extracted by a local diagonal SVD
  (diagnostics only) or deprecated in favour of CTM-only observables.
- A FLOPs/memory cost model for the bond-environment build at the target D.
**Validation:** the revival-rise lag (trajectory RMS) decreases vs Stage-B simple update at
fixed D, monotonically with χ; energy-drift tripwire bounded.
**Risk:** high — research-grade; the indefinite-metric objective is a correctness gap, not
engineering.

### Stage D — full five-site star patch ALS against the CTM environment
**Deliver:** optimize the five updated cores (center + 4 reduced leaves) against the full
plus-star CTM environment (4 external double-layer legs per arm + closed cross interior),
using the exact rank≤2 M7 factorization for the gate. **Cost-gated** against the cheaper
2×2-cluster update (#5, contracts nearest loops exactly, sidesteps env inversion) as a
parallel track.
**Validation:** outperforms Stage C on the trajectory RMS at the entangled peak; passes the
ED-free reverse-echo / energy-drift validators where ED is unavailable (5×5/6×6).
**Risk:** highest. Must NOT collapse the cross to uncoupled bonds (would silently lose the
multi-site correlation the cross creates).

## Validation plan (physics gates, per stage)

| Gate | Stage | Pass/fail criterion | Catches |
|---|---|---|---|
| Native D=2 single-star round-trip vs **dense ncon of the same native tensors**, anisotropic seed (distinct √λ on right vs up) | B | 1-site density & 1-/2-site RDMs match dense to atol 1e-9; FAIL if any density off >1e-8 | transpose/reflection, **leg-routing swap** (invisible at D=1 because the four P↓ legs are identical), √λ-vs-λ mismatch |
| Convention-B equivalence | B | measure via `measurement_peps` vs reabsorb-λ-into-Γ then contract: density/blockade agree to 1e-10 | silent bond-power double-count; future `absorb_s` change |
| Coordinate-reflection trap (single-row-excited, NOT reflection-symmetric) | B | row-resolved density profile matches an explicit `row=Ly−y+1` reference; FAIL if mirrored | `row=Ly−y+1` vs `row=y`, N↔S swap |
| Native vs ED trajectory RMS, 4×4 Néel, [0,2.8], D=2,3,4, **exact oracle, no CTM** | B | rms[D] < 3.5e-2; monotone rms3≤rms2+2e-3, rms4≤rms3+2e-3 (mirror `test_d_convergence.jl`) | native D>1 evolution correctness |
| boundary == dense at **D=5** (and a small-cell D=6) | B/validation | `bnd ≈ dense` atol 1e-9 | the sole D≥5 oracle is currently unchecked at D≥5 |
| CTM residual = real η (not `truncation_error`); χ-sufficiency | C | trust gates on η to tol; environment error ≪ truncation error | "trusted" but unconverged environment driving truncation/ranking |
| Stage-C lag reduction | C | trajectory RMS < Stage-B at fixed D; monotone in χ | the actual mean-field error |
| ED-free reverse-echo R(t) + energy-drift E(t)−E(0) | C/D | bounded drift; reuse `reverse_evolve!`/`exact_pxp_energy_density_finite` | validation at 5×5/6×6 where ED is unavailable |
| Snake-MPO path-independence (only if ever wired in) | C | two snake orderings WITH inter-tensor truncation agree to tol; else reject | path-bias trap |

Cross-cutting: production `validate_pxp_convergence` aggregates a single-time **max**;
the trajectory-RMS methodology lives only in `test_d_convergence.jl` (extended gate, run
by `extended-ctm.yml`). Promote trajectory-RMS into production aggregation.

## Risk assessment

- **Coordinate conventions** — *correct now* (numerically verified) but **fragile**: every
  native gate is D=1, where transpose/reflection/leg-routing-swap all act trivially, and the
  four identical P↓ neighbour legs make a routing swap invisible. Mitigation: the anisotropic
  D=2 and reflection-trap gates above.
- **Gauge/normalization** — two engines use different λ conventions (legacy full-λ, native
  √λ); the converter does the √-split. Compare **observables/CTM summaries, never raw
  tensors** (house rule). Convention-B √λ-in-Γ is unverified numerically at D>1.
- **CTM convergence** — residual mislabeled (truncation_error), gauge not gated, χ=8
  flatters. Fix before any CTM environment drives truncation or energy ranking.
- **Overfitting to D=1** — the dominant systemic risk; the entire suite's native parity is
  gauge-trivial. `verify_m4.jl` (cited as the D>1 CTM-parity proof) is D=1-only.
- **Accidentally changing the model** — keep the five-site cross; never substitute generic
  2-site/`time_evolve`; if the M7 snake is ever wired in, justify it as an *exact*
  factorization and gate path-independence.
- **Unexported PEPSKit internals** (`absorb_s`, `flip_svd`, `truncrank`) — load-error risk
  (safe, not silent); shim + contract test; keep `PEPSKit = "0.7.0"` pinned.
- **Exact oracle at D≥5** — boundary path unverified at the regime where it is the only
  option, and the D=5,6 trajectory is gated on an unimplemented variational
  boundary-MPS recompression.

## Recommended next PR

**Native D>1 validation oracle** (small, self-contained, no new physics; the prerequisite
that makes the existing native backend trustworthy and unblocks everything downstream). See
`docs/superpowers/plans/2026-06-04-native-d-gt-1-validation-oracle.md`. *Not* a rewrite,
*not* CTM-aware evolution — just the gauge-invariant D>1 oracle + the three D>1 traps.

## Verification gates (commands)

```
julia --project=. -e 'using Pkg; Pkg.test()'
SQUAREPXP_EXTENDED_TESTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
```

## Exit criteria

- InfinitePEPS is the only state type; `SUWeight`/λ documented as update-aid + mean-field
  diagnostic; CTMRG is the trusted observable for native states.
- Native D>1 dynamics validated against the exact oracle by the trajectory metric.
- One CTM-aware truncation milestone (Stage C) shows a measurable revival-rise-lag
  reduction at fixed D.
- ITensors/Strided dropped; the dual-engine selector removed.
