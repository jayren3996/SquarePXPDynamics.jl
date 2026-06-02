# Stage-2 code map + diagnosis: bond truncation / gauge / environment (2026-06-02)

Read-only map of the evolution truncation+gauge pipeline, to ground the Stage-2
"good bond truncation + recover good condition with proper regauging" goal.

## What the evolution actually does (it is correct simple-update)
`project_star!` (src/StarSimpleUpdate.jl:453) runs a QR-reduced five-site star
update per center:
1. `_absorb_star_weights` (:250) — each leaf absorbs the FULL λ on its
   center-link **and** its 3 external links (`absorb_link_weight`,
   SquareIPEPS.jl:386 → `T * diag(λ)`).
2. `_qr_reduce_leaves` (:268) — `factorize(...; ortho="left", which_decomp="qr")`
   splits each absorbed leaf into Q (external/environment legs) · R (physical).
3. gate · center · R-factors → reduced core `theta` (:480–487).
4. `_split_reduced_theta` (:307) — sequential SVDs along `split_order`; each:
   - `_svd_with_rel_floor` (:235) applies maxdim/cutoff + the rel_floor relative
     SV condition floor (keep σ ≥ rel_floor·σ_max; caps cond at 1/rel_floor).
   - `_normalized_link_weight_from_singular_values` (:186) → `scale=‖σ‖₂`,
     `λ=σ/scale`; core re-scaled `core = V * scale` (:342) to carry the norm.
5. `_reconstruct_leaf` (:364) — `Q·active`, then `deabsorb_link_weight` the 3
   external λ back out (SquareIPEPS.jl:402 → `T * diag(inv(λ))`, with atol guard).
6. `_commit_star_update!` (:392) — write Γ tensors + new λ, track log-norm.

Weights live in `SquareIPEPSState.link_weights::Dict{BondKey,Vector{Float64}}`
(SquareIPEPS.jl:37), normalized to unit ‖·‖₂, init `[1,0,…,0]`. Hot-path read
`_link_weight_view` (:283).

## Key diagnosis
- Absorbing the **full** λ on environment (external) legs is the textbook
  simple-update mean-field environment (super-orthogonalization: λ² is the
  approximate environment). This is NOT a bug; the evolution is correct.
- Measurement ALREADY does the correct symmetric √λ readout:
  PEPSKitMeasurements.jl:751 splits each bond λ as √λ⊗√λ onto endpoints before
  CTMRG. So the Γ-λ→PEPS conversion for observables is right.
- There is **no regauging code** in evolution (`gauge` field hardcoded `:simple`,
  SquareIPEPS.jl:39); no `regauge`/`canonical`/full-update path exists.
- The only conditioning knob today is `rel_floor` (default 1e-4) — a crude,
  regime-dependent cap that has a bad pocket at the revival.

=> The D=4-worse-than-D=2/3 non-monotonicity AT THE REVIVAL is the inherent
mean-field limitation of simple update on a loopy lattice (chi-sweep already
proved it is an evolution, not a measurement, limit). Closing it requires a
REAL environment inside the evolution loop, not a better knob.

## Stage-2 options, cheapest → most-correct
1. **rel_floor auto / Vidal canonicalization polish** (small): re-canonicalize
   the Γ-λ form between sweeps so λ are true Schmidt values (gauge fix on the
   tree-approx), and/or make rel_floor adaptive. Cheap, may smooth conditioning,
   will NOT remove the mean-field error at the revival.
2. **Loop/cluster-corrected simple update** (medium): include nearest-neighbour
   loop tensors (e.g. simple-update with loop gas / "gauge fixing" à la
   Evenbly, or nth-order cluster) — partial environment, no full CTM.
3. **Full update** (large): use the CTMRG environment (already built for
   measurement) inside the truncation — environment-aware SVD of the bond.
   This is the principled fix for the revival non-monotonicity but is a
   research-grade feature (env construction per bond/sweep, cost, stability).

## Exact hook points for an environment-aware truncation
- `project_star!`:477 — precompute env context before absorption.
- `_qr_reduce_leaves`:276 — where the mean-field (orthogonality) assumption
  enters; swap in an environment-weighted reduction.
- `_split_reduced_theta`:334 — reweight σ by the environment before
  normalize/truncate (the core of full-update truncation).
- `_split_reduced_theta`:342 — adjust `scale` to conserve norm under env weights.
- `_reconstruct_leaf`:375 — match the deabsorption to the chosen gauge.

CTM environment already available via `pepskit_ctmrg_context` /
`to_pepskit_infinitepeps` (PEPSKitMeasurements.jl) — the full-update path can
reuse it rather than building a new contractor.

## DECISIVE EXPERIMENT (2026-06-02): regauging is NOT the Stage-2 fix

Implemented `canonicalize_simple!` (StarSimpleUpdate.jl) — Vidal-gauge
super-orthogonalization via idle (step=0, unprojected, identity-gate) star
sweeps. VERIFIED correct: on an entangled 3x3 state (bond entropy 0.383) it is
state-preserving to machine precision (|Δdensity| 1.4e-13, fidelity exactly 1.0),
idempotent (7 sweeps → converges), and meaningfully regauges (87–100 sweeps/step
during evolution; stored λ shift so bond entropy moves, e.g. 0.3835→0.3633 — the
plain simple-update λ slightly OVER-report entanglement; canonical values are
the true Schmidt entropy).

Measured effect on dynamics (3x3 PXP, :down, dt=0.02, exact_finite vs ED):
- Short time t=0.2, D=4: rel_floor=0 err 9.7e-3 → +regauge 8.7e-3 (only ~10%);
  rel_floor=1e-4 gives 2.5e-5. Regauge does NOT rescue rel_floor=0.
- Entangled t=0.5: rel_floor=1e-4 plain vs +regauge are IDENTICAL to 4–5 sig
  figs (D2 1.2097e-3 both; D3 1.237e-3 vs 1.256e-3; D4 7.33e-4 vs 7.26e-4)
  despite 87–100 regauge sweeps/step. The D=3 wrinkle (f[3]>f[2]) is NOT fixed.
  rel_floor=0+regauge mitigates the D=4 blowup (3.07e-2→2.23e-2) but stays bad.

**Conclusion (physically expected, now empirical):** regauging is a pure gauge
transformation — it leaves the represented wavefunction unchanged, so it cannot
change what the next Trotter step truncates. The simple-update truncation quality
is governed by the mean-field (λ²) ENVIRONMENT, which is identical in every gauge.
Therefore proper regauging — the literal Stage-2 "recover good condition" ask and
the roadmap's hypothesized fix — is EMPIRICALLY RULED OUT as sufficient. The
remaining real fix is environment-aware truncation (full update with the CTM bond
environment). `rel_floor` stays the practical conditioning knob until then.

`canonicalize_simple!` is kept as a correct, tested primitive (accurate Schmidt /
entanglement diagnostics, and the canonicalization a full-update step will need),
but is NOT wired into `evolve!` — it is not a dynamics fix.
