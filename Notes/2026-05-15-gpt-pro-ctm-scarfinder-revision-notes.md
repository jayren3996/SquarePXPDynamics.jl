# GPT-Pro CTM ScarFinder Revision Notes

Date saved: 2026-05-15

This note captures the review takeaway and a practical revision sketch for hardening the current square-lattice PXP/iPEPS code before using ScarFinder scores for energy-targeted or physics-quality ranking.

## Executive Takeaway

The repo is in a healthy prototype state for continuing S6b-style ScarFinder candidate ranking and logging, but it is not ready for energy-targeted ScarFinder decisions.

The architecture is mostly sound:

- Dense square-star PXP remains the source of truth.
- Five-color scheduling is deterministic.
- `project_star!` is the low-level one-star update primitive.
- `evolve!` orchestrates Trotter layers over `project_star!`.
- `scarfinder!` stays above evolution and measurement.
- Simple/local observables are documented as diagnostics, not CTMRG-quality physics.

The main risk is numerical trust, not broad code organization.

## Main Blockers

1. CTMRG diagnostics are not production-like yet.
   - Raw PEPSKit `info` is preserved, but not interpreted.
   - Public summaries do not expose convergence status, iterations, residuals, chi, tolerance, or acceptance policy.
   - There is no context/state compatibility guard.
   - There is no finite-chi sensitivity workflow.

2. D=2 algorithm tests are still mostly smoke tests.
   - D=1 dense-reference coverage is good.
   - D=2 needs gauge-invariant regression tests, not only finite/norm checks.

3. Repeated-update norm and gauge control need hardening.
   - The QR/SVD path normalizes each new lambda spectrum and pushes scale into the remaining core.
   - This reconstructs a single update, but long ScarFinder runs can accumulate uncontrolled Gamma tensor norm factors.

4. Small validation/API gaps remain.
   - `PEPSKitCTMRGParams` should reject nonfinite `tol`.
   - `square_star_basis_allowed(bits)` should validate that user-supplied bits are only `0` or `1`.
   - `project_star!(psi; maxdim > psi.maxdim)` needs an explicit policy.

## Recommended Milestone

Choose C: improve CTMRG diagnostics and convergence, and make CTM measurements more production-like.

Candidate ranking and logging can continue as diagnostic plumbing, but energy targeting should wait until CTM energy and CTM convergence reporting are harder to misuse.

## Revision Plan

### Phase 0: Baseline and Policy Decisions

Run the current test suite before changing behavior:

```bash
julia --project -e 'using Pkg; Pkg.test()'
SQUAREPXP_EXTENDED_TESTS=1 julia --project -e 'using Pkg; Pkg.test()'
```

Inspect the actual PEPSKit `info` structure returned by `leading_boundary`.

Decide these policies before implementation:

- Stale CTM context behavior: throw, warn, or document as caller error. Recommended: throw.
- `project_star!(maxdim > psi.maxdim)` behavior: either allow growth and update/clarify metadata, or reject it. Recommended: allow growth only if `psi.maxdim` is documented as a default cap rather than an invariant.
- Unconverged CTM ranking behavior: either flag candidates as CTM-untrusted or require explicit opt-in to rank with unconverged CTM values.

### Phase 1: Low-Risk Validation Fixes

Likely files:

- `src/SquarePXP.jl`
- `src/PEPSKitMeasurements.jl`
- `test/test_square_pxp.jl`
- `test/test_pepskit_measurements.jl`
- `test/test_star_simple_update.jl`

Add failing tests first:

- `square_star_basis_allowed((0, 2, 1, 1, 1))` throws `ArgumentError`.
- Similar invalid-bit cases throw.
- `PEPSKitCTMRGParams(4, Inf, 10, 0)` throws.
- `PEPSKitCTMRGParams(4, NaN, 10, 0)` throws.
- Chosen `maxdim` growth/staleness policy is enforced.

Then implement the minimal validation changes.

### Phase 2: CTMRG Diagnostics and Context Safety

Likely files:

- `src/PEPSKitMeasurements.jl`
- `src/ScarFinder.jl`
- `test/test_pepskit_measurements.jl`
- `test/test_scarfinder.jl`

Add a structured diagnostics type, for example:

```julia
struct CTMRGDiagnostics
    chi::Int
    tol::Float64
    maxiter::Int
    iterations::Union{Int,Nothing}
    residual::Union{Float64,Nothing}
    converged::Union{Bool,Nothing}
    accepted::Bool
end
```

The exact fields should match what PEPSKit exposes through `info`, but the public API should at least report:

- chi,
- tolerance,
- max iterations,
- actual iterations if available,
- residual/error if available,
- convergence status if available,
- whether the measurement was accepted for downstream ranking.

Add a context/state freshness guard:

- Store a fingerprint or version in `PEPSKitMeasurementContext`.
- Check it in CTM measurement calls.
- Throw if a context created from one state is used with another state or with a mutated state.

Possible implementation options:

- Add a mutation/version counter to `SquareIPEPSState`, incremented by `project_star!`.
- Or store a lightweight fingerprint derived from unit-cell shape, tensor/link index structure, and link-weight dimensions.

The version-counter approach is cleaner if mutable state changes are centralized.

### Phase 3: ScarFinder CTM Logging and Ranking Hardening

Likely files:

- `src/ScarFinder.jl`
- `test/test_scarfinder.jl`

Add CTM diagnostic fields to candidate summaries and logs.

Expected behavior:

- CTM callback summaries include convergence metadata.
- ScarFinder logs preserve CTM trust status.
- Energy-oriented ranking refuses or flags unconverged CTM diagnostics.
- Existing simple diagnostic ranking remains available, but is documented as a diagnostic sorting key.

Do not introduce true energy targeting in this phase. First make current ranking honest.

### Phase 4: D=2 Regression Tests

Likely files:

- `test/test_star_simple_update.jl`
- `test/test_observables_evolved.jl`

Add deterministic D=2 fixtures:

- Seeded random complex Gamma tensors.
- Positive normalized lambda vectors.
- A small measurement panel using gauge-invariant/simple observables:
  - selected-site density,
  - nearest-neighbor density,
  - selected star expectations.

Tests to add:

- Zero-step `projected=false` identity on a nontrivial D=2 state.
- Split-order equivalence using observable comparisons, not raw tensor equality.
- If cheap, an exactly representable update regression comparing local star expectations before and after.

These tests should catch leg-ordering, lambda absorption/deabsorption, and SVD replacement bugs that D=1 cannot catch.

### Phase 5: Norm-Scale Tracking

Likely files:

- `src/SquareIPEPS.jl`
- `src/StarSimpleUpdate.jl`
- `src/IPEPSEvolution.jl`
- `test/test_star_simple_update.jl`

Add one of:

- a `log_norm::Float64` or similar normalization ledger on the state, updated from split `norm_factors`; or
- explicit Gamma tensor rescaling with recorded scale.

Recommended first approach: add a ledger. It is less invasive than rescaling tensors and easier to test.

Tests:

- Run 100-500 small D=2 updates.
- Assert Gamma norms remain bounded if rescaling is used, or that log norm tracks accumulated scaling if a ledger is used.
- Check that gauge-invariant observables are unchanged by pure rescaling.

### Phase 6: CTM Energy Behavior Tests

Likely file:

- `test/test_pepskit_measurements.jl`

Replace source-inspection tests with behavior tests:

- Compare `pxp_energy_density_ctm(psi, ctx)` with the average of `star_expectation_ctm(psi, c, Hstar, ctx)` over all centers using the same context.
- Run at least on a product or D=1 state.
- Put a short-evolved-state version behind `SQUAREPXP_EXTENDED_TESTS=1` if runtime is high.

This protects the true CTM-backed energy path without brittle regex/source checks.

Status:

- Added default product-state coverage comparing `pxp_energy_density_ctm(psi, ctx)`
  against the average of per-center `star_expectation_ctm(psi, c, Hstar, ctx)`.
- Added extended short-evolved D=1 coverage behind `SQUAREPXP_EXTENDED_TESTS=1`
  using one PEPSKit CTMRG context and modest `PEPSKitCTMRGParams(2, 1e-6, 20, 0)`.
- Deferred a D=2 CTM energy stress test because the extended PEPSKit CTMRG
  measurement file already takes about five minutes on the current setup.

### Phase 7: Documentation

Likely file:

- `README.md`

Add a small CTM example:

```julia
params_ctm = PEPSKitCTMRGParams(8, 1e-8, 100, 0)
ctx = pepskit_ctmrg_context(psi; params=params_ctm)
energy = pxp_energy_density_ctm(psi, ctx)
```

Immediately state:

- Users must inspect CTMRG diagnostics before using CTM values for ranking.
- Users should repeat CTM measurements at multiple chi values before trusting energy comparisons.
- A `PEPSKitMeasurementContext` belongs to the exact state used at creation.
- If `psi` is mutated by `evolve!` or `project_star!`, the old context is stale.
- `ScarFinderCandidateScore.score` is a diagnostic sorting key, not a physics-quality energy target.

## Suggested First PR Slice

Start with a small revision containing:

1. Validation fixes.
2. CTMRG diagnostics type.
3. Context freshness check.
4. CTM metadata in summaries/logs.
5. Tests for all of the above.

This gives immediate safety without mixing in D=2 algorithm regressions or norm-ledger design.

Then follow with:

1. D=2 gauge-invariant regression tests.
2. Norm-scale tracking.
3. CTM energy behavior tests.
4. README updates.

## Deferred Work

Defer the following until CTM measurements and D=2 update regressions are trustworthy:

- full-update gauge fixing,
- true energy targeting/correction,
- production ScarFinder validation,
- PEPSKit/TensorKit package-extension split,
- threaded disjoint-star updates,
- array-backed hot-path kernels,
- broad module splitting.
