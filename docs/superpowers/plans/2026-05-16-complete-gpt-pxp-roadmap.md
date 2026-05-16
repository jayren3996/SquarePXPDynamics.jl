# Complete GPT PXP Roadmap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Finish the remaining GPT review items by making ScarFinder consume trusted CTM-backed physics objectives, adding the missing observables and convergence reports, making runs reproducible and persistent, and documenting remaining CTM-aware update boundaries.

**Architecture:** Keep the existing custom five-site simple-update path as the mutation engine, and promote `PXPValidation.jl` plus `PEPSKitMeasurements.jl` into the trusted measurement layer consumed by `ScarFinder.jl`. Work in slices that each leave the package passing `Pkg.test()`: first backend/objective abstractions, then observables, reproducibility, convergence/error budgets, persistence, reverse-evolution tests, projection semantics, and finally CTM-aware update readiness APIs.

**Tech Stack:** Julia 1.12, ITensors, PEPSKit, TensorKit, EDKit, JSON3, existing package tests under `test/runtests.jl`.

---

## Scope And File Map

- Modify `src/ScarFinder.jl`: add measurement backend dispatch, objective types, trusted CTM ranking, candidate metadata, and persistence hooks without moving low-level tensor logic into ScarFinder.
- Modify `src/PXPValidation.jl`: add convergence/error-budget report types, reverse-evolution validation helpers, and machine-readable report writers.
- Modify `src/PEPSKitMeasurements.jl`: add deterministic CTMRG initialization controls and CTM-backed observables needed by ScarFinder.
- Modify `src/Observables.jl`: add simple analogs of imbalance and structure-factor-like checkerboard diagnostics for tests and fallback logs.
- Modify `src/SquarePXP.jl`: clarify and optionally expose `P * U * P` projection semantics while preserving current `P * U` default behavior.
- Modify `src/SquarePXPDynamics.jl`: re-export every new public type/function and satisfy public-doc tests.
- Create `src/ScarFinderStores.jl` only if persistence grows beyond a small writer in `ScarFinder.jl`; otherwise keep persistence local to `ScarFinder.jl`.
- Modify tests: `test/test_scarfinder.jl`, `test/test_pxp_validation.jl`, `test/test_pepskit_measurements.jl`, `test/test_observables.jl`, `test/test_square_pxp.jl`, `test/test_public_docs.jl` indirectly through docstrings.
- Modify docs: `README.md`, `memory/mid_term/decision_log.md`.

## Task 1: ScarFinder Measurement Backends

**Files:**
- Modify: `src/ScarFinder.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_scarfinder.jl`

- [ ] **Step 1: Write backend construction tests**

Add this testset near the existing CTM callback tests in `test/test_scarfinder.jl`:

```julia
@testset "ScarFinder measurement backends construct and validate" begin
    simple = SimpleBackend()
    @test simple isa MeasurementBackend

    params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    policy = CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4)
    trusted = TrustedCTMBackend(params, policy)

    @test trusted isa MeasurementBackend
    @test trusted.params === params
    @test trusted.policy === policy
    @test_throws ArgumentError TrustedCTMBackend((), policy)
end
```

- [ ] **Step 2: Run test and confirm missing symbols**

Run:

```bash
julia --project=. test/runtests.jl test_scarfinder.jl
```

Expected: FAIL with `UndefVarError: SimpleBackend not defined`.

- [ ] **Step 3: Add backend types and measurement helper**

In `src/ScarFinder.jl`, extend imports:

```julia
using ..PXPValidation: TrustedCTMMeasurement, measure_ctm_trusted
using ..CTMTrust: CTMTrustPolicy
using ..PEPSKitMeasurements: PEPSKitCTMRGParams, measure_ctm
```

Add exports:

```julia
export MeasurementBackend, SimpleBackend, TrustedCTMBackend, measure_scarfinder
```

Add these definitions after `ScarFinderParams`:

```julia
abstract type MeasurementBackend end

"""Development backend using local/simple-update observables only."""
struct SimpleBackend <: MeasurementBackend end

"""Trusted finite-chi CTM backend for production ScarFinder diagnostics."""
struct TrustedCTMBackend <: MeasurementBackend
    params::Tuple
    policy::CTMTrustPolicy
    measure::Any

    function TrustedCTMBackend(
        params,
        policy::CTMTrustPolicy = CTMTrustPolicy();
        measure = measure_ctm,
    )
        collected = Tuple(params)
        isempty(collected) &&
            throw(ArgumentError("TrustedCTMBackend requires at least one CTMRG parameter point"))
        all(p -> p isa PEPSKitCTMRGParams, collected) ||
            throw(ArgumentError("TrustedCTMBackend params must be PEPSKitCTMRGParams"))
        return new(collected, policy, measure)
    end
end

"""Measure `psi` with a ScarFinder measurement backend."""
measure_scarfinder(psi::SquareIPEPSState, ::SimpleBackend) = measure_simple(psi)

function measure_scarfinder(psi::SquareIPEPSState, backend::TrustedCTMBackend)
    return measure_ctm_trusted(
        psi;
        params = backend.params,
        policy = backend.policy,
        measure = backend.measure,
    )
end
```

- [ ] **Step 4: Export from the top-level module**

In `src/SquarePXPDynamics.jl`, import and export:

```julia
using .ScarFinder:
    MeasurementBackend,
    SimpleBackend,
    TrustedCTMBackend,
    measure_scarfinder,
    ScarFinderParams,
    ScarFinderCandidateScore,
    ScarFinderIteration,
    ScarFinderResult,
    rank_scarfinder_candidates,
    write_scarfinder_log,
    scarfinder!

export MeasurementBackend, SimpleBackend, TrustedCTMBackend, measure_scarfinder
```

- [ ] **Step 5: Run focused tests**

Run:

```bash
julia --project=. test/runtests.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ScarFinder.jl src/SquarePXPDynamics.jl test/test_scarfinder.jl
git commit -m "feat: add ScarFinder measurement backends"
```

## Task 2: Physics Objectives And Trusted Ranking

**Files:**
- Modify: `src/ScarFinder.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_scarfinder.jl`

- [ ] **Step 1: Write objective tests**

Add this testset:

```julia
@testset "ScarFinder physics objectives score candidates" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    objective = CompositeObjective(;
        revival = RevivalObjective(:sublattice_imbalance, 1.0),
        blockade_weight = 10.0,
        truncation_weight = 2.0,
        finite_chi_weight = 3.0,
        entropy_weight = 0.5,
    )

    result = scarfinder!(psi, params; objective)
    score = only(rank_scarfinder_candidates(result; diagnostics = :simple))

    @test score.objective_name == "CompositeObjective"
    @test isfinite(score.score)
    @test score.revival_strength !== nothing
    @test score.finite_chi_drift === nothing
end
```

- [ ] **Step 2: Run test and confirm missing objective symbols**

Run:

```bash
julia --project=. test/runtests.jl test_scarfinder.jl
```

Expected: FAIL with `UndefVarError: CompositeObjective not defined`.

- [ ] **Step 3: Add objective types**

In `src/ScarFinder.jl`, export:

```julia
export ScarFinderObjective, RevivalObjective, TargetEnergyObjective
export LowVarianceObjective, CompositeObjective
```

Add:

```julia
abstract type ScarFinderObjective end

struct RevivalObjective <: ScarFinderObjective
    observable::Symbol
    weight::Float64
    function RevivalObjective(observable::Symbol = :sublattice_imbalance, weight::Real = 1.0)
        observable in (:sublattice_imbalance, :density) ||
            throw(ArgumentError("revival observable must be :sublattice_imbalance or :density"))
        w = Float64(weight)
        isfinite(w) && w >= 0 || throw(ArgumentError("revival weight must be finite and nonnegative"))
        return new(observable, w)
    end
end

struct TargetEnergyObjective <: ScarFinderObjective
    target::Float64
    weight::Float64
end

struct LowVarianceObjective <: ScarFinderObjective
    weight::Float64
end

struct CompositeObjective <: ScarFinderObjective
    revival::Union{Nothing,RevivalObjective}
    target_energy::Union{Nothing,TargetEnergyObjective}
    low_variance::Union{Nothing,LowVarianceObjective}
    blockade_weight::Float64
    truncation_weight::Float64
    finite_chi_weight::Float64
    entropy_weight::Float64

    function CompositeObjective(;
        revival::Union{Nothing,RevivalObjective} = RevivalObjective(),
        target_energy::Union{Nothing,TargetEnergyObjective} = nothing,
        low_variance::Union{Nothing,LowVarianceObjective} = nothing,
        blockade_weight::Real = 100.0,
        truncation_weight::Real = 10.0,
        finite_chi_weight::Real = 10.0,
        entropy_weight::Real = 1.0,
    )
        weights = Float64.((blockade_weight, truncation_weight, finite_chi_weight, entropy_weight))
        all(w -> isfinite(w) && w >= 0, weights) ||
            throw(ArgumentError("objective weights must be finite and nonnegative"))
        return new(revival, target_energy, low_variance, weights...)
    end
end
```

- [ ] **Step 4: Extend score records**

Add fields to `ScarFinderCandidateScore` after `score::Float64`:

```julia
objective_name::String
revival_strength::Union{Nothing,Float64}
finite_chi_drift::Union{Nothing,Float64}
energy_variance_proxy::Union{Nothing,Float64}
```

Update every `ScarFinderCandidateScore(...)` constructor call with:

```julia
score,
String(nameof(typeof(objective))),
revival_strength,
finite_chi_drift,
nothing,
```

Use a local `objective::ScarFinderObjective = CompositeObjective()` argument in `_candidate_score`.

- [ ] **Step 5: Implement scoring**

Add helpers:

```julia
function _imbalance(obs)
    return obs.density_even - obs.density_odd
end

function _revival_strength(obs, objective::RevivalObjective)
    objective.observable === :sublattice_imbalance && return abs(_imbalance(obs))
    objective.observable === :density && return abs(obs.density)
    throw(ArgumentError("unsupported revival observable"))
end

_finite_chi_drift(::Nothing) = nothing
function _finite_chi_drift(trusted::TrustedCTMMeasurement)
    vals = (
        trusted.trust.finite_chi_density_delta,
        trusted.trust.finite_chi_blockade_delta,
        trusted.trust.finite_chi_energy_delta,
    )
    present = collect(skipmissing(v === nothing ? missing : v for v in vals))
    return isempty(present) ? nothing : maximum(abs, present)
end

function _score_value(obs, log::EvolutionLog, mean_bond_entropy, max_bond_entropy, objective::CompositeObjective; trusted_ctm = nothing)
    revival = objective.revival === nothing ? nothing : _revival_strength(obs, objective.revival)
    finite_chi = _finite_chi_drift(trusted_ctm)
    entropy_penalty = max_bond_entropy === nothing ? 0.0 : max_bond_entropy
    score = obs.blockade_violation * objective.blockade_weight +
            log.max_truncerr * objective.truncation_weight +
            entropy_penalty * objective.entropy_weight
    finite_chi === nothing || (score += finite_chi * objective.finite_chi_weight)
    revival === nothing || (score -= objective.revival.weight * revival)
    objective.target_energy === nothing ||
        (score += objective.target_energy.weight * abs(obs.pxp_energy_density - objective.target_energy.target))
    return score, revival, finite_chi
end
```

- [ ] **Step 6: Thread objective through `scarfinder!`**

Add keyword:

```julia
objective::ScarFinderObjective = CompositeObjective()
```

to both `scarfinder!` methods. Pass it into `_candidate_score`.

- [ ] **Step 7: Update log writers**

Add CSV/JSON fields:

```julia
"objective_name", "revival_strength", "finite_chi_drift", "energy_variance_proxy"
```

and include corresponding values in `_write_csv_log` and `_score_json`.

- [ ] **Step 8: Run focused tests**

```bash
julia --project=. test/runtests.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add src/ScarFinder.jl src/SquarePXPDynamics.jl test/test_scarfinder.jl
git commit -m "feat: add ScarFinder physics objectives"
```

## Task 3: Trusted CTM Backend As Production ScarFinder Path

**Files:**
- Modify: `src/ScarFinder.jl`
- Test: `test/test_scarfinder.jl`

- [ ] **Step 1: Write trusted-backend acceptance test**

```julia
@testset "ScarFinder trusted CTM backend gates ranking" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    ctm_params = (
        PEPSKitCTMRGParams(2, 1e-5, 4, 0),
        PEPSKitCTMRGParams(4, 1e-6, 4, 0),
    )
    backend = TrustedCTMBackend(
        ctm_params,
        CTMTrustPolicy(2, true, 1e-2, 1e-3, 1e-2, 1e-4);
        measure = (state; params) -> CTMObservableSummary(
            0.1 + params.chi * 1e-5,
            0.12,
            0.08,
            params.chi * 1e-6,
            -0.01,
            CTMRGDiagnostics(params.chi, params.tol, params.maxiter, 3, params.tol / 10, true, true),
        ),
    )

    result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        measurement = backend,
        ctm_every = 1,
        require_trusted_ctm = true,
    )

    ranked = rank_scarfinder_candidates(result; diagnostics = :ctm, require_ctm_trusted = true)
    @test length(ranked) == 1
    @test ranked[1].ctm_trusted === true
    @test ranked[1].ctm_trust_reason == "trusted"
    @test ranked[1].finite_chi_drift !== nothing
end
```

- [ ] **Step 2: Add score trust fields**

Add to `ScarFinderCandidateScore`:

```julia
ctm_trusted::Union{Nothing,Bool}
ctm_trust_reason::Union{Nothing,String}
```

Populate `nothing` for simple scores, and populate from `TrustedCTMMeasurement.trust` for trusted CTM scores.

- [ ] **Step 3: Accept trusted CTM measurement results**

Add `_candidate_score` method:

```julia
function _candidate_score(
    iteration::Int,
    diagnostics::Symbol,
    accepted::Bool,
    reject_reason::Union{Nothing,String},
    log::EvolutionLog,
    trusted::TrustedCTMMeasurement,
    objective::ScarFinderObjective,
    correction_accepted::Union{Nothing,Bool} = nothing,
    correction_energy_before::Union{Nothing,Float64} = nothing,
    correction_energy_after::Union{Nothing,Float64} = nothing,
)
    score = _candidate_score(
        iteration,
        diagnostics,
        accepted,
        reject_reason,
        log,
        trusted.measurement,
        objective,
        correction_accepted,
        correction_energy_before,
        correction_energy_after;
        trusted_ctm = trusted,
    )
    return score
end
```

- [ ] **Step 4: Replace callback-only measurement flow**

Add keyword to `scarfinder!`:

```julia
measurement::MeasurementBackend = SimpleBackend()
require_trusted_ctm::Bool = false
```

Inside scheduled CTM block, choose:

```julia
ctm_obs = ctm_callback === nothing ? measure_scarfinder(psi, measurement) : ctm_callback(psi, n, simple_score)
```

If `require_trusted_ctm` and `ctm_obs isa TrustedCTMMeasurement` and `!ctm_obs.trust.trusted`, set `accepted = false` and `reason = "trusted CTM policy rejected iteration"` before building scores.

- [ ] **Step 5: Extend ranking filter**

Add keyword:

```julia
require_ctm_trusted::Bool = false
```

with validation `require_ctm_trusted && diagnostics !== :ctm && throw(...)`, and filter with `score.ctm_trusted === true`.

- [ ] **Step 6: Run tests**

```bash
julia --project=. test/runtests.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add src/ScarFinder.jl test/test_scarfinder.jl
git commit -m "feat: make ScarFinder consume trusted CTM measurements"
```

## Task 4: CTM-Backed Scar Observables

**Files:**
- Modify: `src/Observables.jl`
- Modify: `src/PEPSKitMeasurements.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_observables.jl`
- Test: `test/test_pepskit_measurements.jl`

- [ ] **Step 1: Add simple observable tests**

```julia
@testset "PXP scar observables simple product limits" begin
    cell = PeriodicSquareUnitCell(10, 10)
    down = product_square_ipeps(cell; state = :down, maxdim = 1)
    checker = checkerboard_square_ipeps(cell; maxdim = 1)

    @test sublattice_imbalance_simple(down) ≈ 0.0 atol = 1e-12
    @test abs(sublattice_imbalance_simple(checker)) ≈ 1.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(down) ≈ 0.0 atol = 1e-12
    @test checkerboard_structure_factor_simple(checker) ≈ 1.0 atol = 1e-12
end
```

- [ ] **Step 2: Implement simple helpers**

In `src/Observables.jl`:

```julia
export sublattice_imbalance_simple, checkerboard_structure_factor_simple

"""Return even-minus-odd Rydberg density from simple/local measurements."""
function sublattice_imbalance_simple(psi::SquareIPEPSState)::Float64
    even, odd = sublattice_densities(psi)
    return even - odd
end

"""Return the squared checkerboard density contrast from simple/local measurements."""
function checkerboard_structure_factor_simple(psi::SquareIPEPSState)::Float64
    return sublattice_imbalance_simple(psi)^2
end
```

- [ ] **Step 3: Add CTM summary fields**

Extend `CTMObservableSummary` with:

```julia
sublattice_imbalance::Float64
checkerboard_structure_factor::Float64
```

Update constructors so legacy five-argument calls compute:

```julia
imbalance = density_even - density_odd
structure = imbalance^2
```

- [ ] **Step 4: Include fields in finite checks and JSON/report data**

Update `_assert_finite_ctm_summary`, `_ctm_summary_data` in `PXPValidation.jl`, ScarFinder score construction, and CTM validation point comparisons if needed. Keep existing density deltas unchanged.

- [ ] **Step 5: Run focused tests**

```bash
julia --project=. test/runtests.jl test_observables.jl test_pepskit_measurements.jl test_pxp_validation.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/Observables.jl src/PEPSKitMeasurements.jl src/PXPValidation.jl src/ScarFinder.jl src/SquarePXPDynamics.jl test/test_observables.jl test/test_pepskit_measurements.jl test/test_pxp_validation.jl test/test_scarfinder.jl
git commit -m "feat: add scar-oriented PXP observables"
```

## Task 5: Deterministic CTMRG Initialization

**Files:**
- Modify: `src/PEPSKitMeasurements.jl`
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_pepskit_measurements.jl`
- Test: `test/test_pxp_validation.jl`

- [ ] **Step 1: Write parameter tests**

```julia
@testset "PEPSKit CTMRG params carry reproducibility seed" begin
    params = PEPSKitCTMRGParams(4, 1e-8, 10, 0; seed = 1234)
    @test params.seed == 1234
    @test PEPSKitCTMRGParams(4, 1e-8, 10, 0).seed === nothing
    @test_throws ArgumentError PEPSKitCTMRGParams(4, 1e-8, 10, 0; seed = -1)
end
```

- [ ] **Step 2: Add `Random` and seed field**

In `src/PEPSKitMeasurements.jl`:

```julia
using Random
```

Extend `PEPSKitCTMRGParams`:

```julia
seed::Union{Nothing,Int}
```

Change constructor signature:

```julia
function PEPSKitCTMRGParams(chi::Integer, tol::Real, maxiter::Integer, verbosity::Integer; seed::Union{Nothing,Integer} = nothing)
```

Validate:

```julia
seed_value = seed === nothing ? nothing : Int(seed)
seed_value === nothing || seed_value >= 0 || throw(ArgumentError("seed must be nonnegative"))
return new(Int(chi), Float64(tol), Int(maxiter), Int(verbosity), seed_value)
```

- [ ] **Step 3: Use seeded CTMRG environment initialization**

In `pepskit_ctmrg_context`:

```julia
rng = params.seed === nothing ? Random.default_rng() : Random.MersenneTwister(params.seed)
env0 = PEPSKit.CTMRGEnv((args...) -> randn(rng, args...), ComplexF64, peps, chi_space)
```

If PEPSKit rejects that callable shape, replace with a private helper:

```julia
_ctm_initializer(::Nothing) = randn
function _ctm_initializer(seed::Int)
    rng = Random.MersenneTwister(seed)
    return (T, dims...) -> randn(rng, T, dims...)
end
```

and call `PEPSKit.CTMRGEnv(_ctm_initializer(params.seed), ComplexF64, peps, chi_space)`.

- [ ] **Step 4: Serialize seeds**

Add `seed = params.seed` to CTM point JSON in `src/PXPValidation.jl` and to CSV output in `write_ctm_validation_csv` if present.

- [ ] **Step 5: Run tests**

```bash
julia --project=. test/runtests.jl test_pepskit_measurements.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/PEPSKitMeasurements.jl src/PXPValidation.jl test/test_pepskit_measurements.jl test/test_pxp_validation.jl
git commit -m "feat: add reproducible CTMRG initialization"
```

## Task 6: Error-Budget And Convergence Reports

**Files:**
- Modify: `src/PXPValidation.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Create: `scripts/pxp_convergence_report.jl`
- Test: `test/test_pxp_validation.jl`

- [ ] **Step 1: Write convergence config tests**

```julia
@testset "PXP convergence report aggregates validation grid" begin
    base = PXPValidationConfig(3; total_time = 0.01, dt = 0.01, measure_every = 1)
    sweep = PXPConvergenceConfig(
        base;
        dt_values = [0.01, 0.005],
        D_values = [1],
        chi_values = Int[],
        cutoff_values = [1e-12],
    )

    report = validate_pxp_convergence(
        sweep;
        ctm_measure = (state; params) -> CTMObservableSummary(0.0, 0.0, 0.0, 0.0, 0.0),
    )

    @test length(report.runs) == 2
    @test report.runs[1].config.dt == 0.01
    @test report.runs[2].config.dt == 0.005
    @test isfinite(report.max_abs_density_error_simple)
end
```

- [ ] **Step 2: Add report types**

In `src/PXPValidation.jl`:

```julia
export PXPConvergenceConfig, PXPConvergenceReport, validate_pxp_convergence

struct PXPConvergenceConfig
    base::PXPValidationConfig
    dt_values::Vector{Float64}
    D_values::Vector{Int}
    chi_values::Vector{Int}
    cutoff_values::Vector{Float64}
end

struct PXPConvergenceReport
    config::PXPConvergenceConfig
    runs::Vector{PXPValidationReport}
    max_abs_density_error_simple::Float64
    max_abs_density_error_ctm::Union{Nothing,Float64}
    all_ctm_trusted::Union{Nothing,Bool}
end
```

- [ ] **Step 3: Implement grid runner**

```julia
function _copy_config(base::PXPValidationConfig; dt = base.dt, maxdim = base.maxdim, cutoff = base.cutoff)
    return PXPValidationConfig(
        base.n;
        total_time = base.total_time,
        dt,
        measure_every = max(1, round(Int, base.measure_every * base.dt / dt)),
        initial_state = base.initial_state,
        point_group = base.point_group,
        use_sparse = base.use_sparse,
        ed_tol = base.ed_tol,
        ed_m_init = base.ed_m_init,
        ed_m_max = base.ed_m_max,
        ed_extend_step = base.ed_extend_step,
        order = base.order,
        maxdim,
        cutoff,
        schedule = base.schedule,
    )
end

function validate_pxp_convergence(config::PXPConvergenceConfig; ctm_measure = measure_ctm)
    runs = PXPValidationReport[]
    for dt in config.dt_values, D in config.D_values, cutoff in config.cutoff_values
        run_config = _copy_config(config.base; dt, maxdim = D, cutoff)
        push!(runs, validate_pxp_ed_ipeps(run_config; ctm_params = nothing, ctm_measure))
    end
    simple_errors = [abs(c.density_error_simple) for r in runs for c in r.comparisons]
    return PXPConvergenceReport(config, runs, maximum(simple_errors), nothing, nothing)
end
```

Add CTM `chi_values` support in the same function once Task 5 seed serialization is complete:

```julia
ctm_params = isempty(config.chi_values) ? nothing :
    Tuple(PEPSKitCTMRGParams(chi, 1e-8, 100, 0) for chi in config.chi_values)
```

- [ ] **Step 4: Add script**

Create `scripts/pxp_convergence_report.jl`:

```julia
using JSON3
using SquarePXPDynamics

base = PXPValidationConfig(
    parse(Int, get(ENV, "SQUAREPXP_CONVERGENCE_N", "3"));
    total_time = parse(Float64, get(ENV, "SQUAREPXP_CONVERGENCE_TOTAL_TIME", "0.02")),
    dt = parse(Float64, get(ENV, "SQUAREPXP_CONVERGENCE_BASE_DT", "0.01")),
)
sweep = PXPConvergenceConfig(
    base;
    dt_values = parse.(Float64, split(get(ENV, "SQUAREPXP_CONVERGENCE_DT", "0.01,0.005"), ",")),
    D_values = parse.(Int, split(get(ENV, "SQUAREPXP_CONVERGENCE_D", "1"), ",")),
    chi_values = Int[],
    cutoff_values = parse.(Float64, split(get(ENV, "SQUAREPXP_CONVERGENCE_CUTOFF", "1e-12"), ",")),
)
report = validate_pxp_convergence(sweep)
out = get(ENV, "SQUAREPXP_CONVERGENCE_OUT", "pxp_convergence_report.json")
write_pxp_convergence_json(report, out)
println(out)
```

- [ ] **Step 5: Run tests and script smoke**

```bash
julia --project=. test/runtests.jl test_pxp_validation.jl test_public_docs.jl
julia --project=. scripts/pxp_convergence_report.jl
```

Expected: tests PASS and script prints `pxp_convergence_report.json`.

- [ ] **Step 6: Commit**

```bash
git add src/PXPValidation.jl src/SquarePXPDynamics.jl scripts/pxp_convergence_report.jl test/test_pxp_validation.jl
git commit -m "feat: add PXP convergence reports"
```

## Task 7: Candidate Persistence

**Files:**
- Modify: `src/ScarFinder.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_scarfinder.jl`

- [ ] **Step 1: Write persistence test**

```julia
@testset "ScarFinder candidate store writes auditable metadata" begin
    cell = PeriodicSquareUnitCell(10, 10)
    trotter = TrotterParams(0.01, 1, :real, true, 1, 1e-12)
    params = ScarFinderParams(0.0, trotter, 1, Inf, Inf, Inf, false)
    dir = mktempdir()
    store = JSONCandidateStore(dir)

    result = scarfinder!(
        product_square_ipeps(cell; state = :down, maxdim = 1),
        params;
        candidate_store = store,
    )

    files = readdir(dir)
    @test "candidate_000001.json" in files
    data = JSON3.read(read(joinpath(dir, "candidate_000001.json"), String))
    @test data.iteration == 1
    @test data.accepted == true
    @test data.state_version isa Integer
    @test data.log_norm isa Number
    @test data.score.diagnostics == "simple"
end
```

- [ ] **Step 2: Add store types**

In `src/ScarFinder.jl`:

```julia
using JSON3

export CandidateStore, NoCandidateStore, JSONCandidateStore

abstract type CandidateStore end
struct NoCandidateStore <: CandidateStore end

struct JSONCandidateStore <: CandidateStore
    directory::String
    function JSONCandidateStore(directory::AbstractString)
        mkpath(directory)
        return new(String(directory))
    end
end
```

- [ ] **Step 3: Add writer**

```julia
store_candidate!(::NoCandidateStore, psi::SquareIPEPSState, iteration::ScarFinderIteration) = nothing

function store_candidate!(store::JSONCandidateStore, psi::SquareIPEPSState, iteration::ScarFinderIteration)
    path = joinpath(store.directory, "candidate_$(lpad(iteration.iteration, 6, '0')).json")
    payload = (;
        iteration = iteration.iteration,
        accepted = iteration.accepted,
        reject_reason = iteration.reject_reason,
        state_version = state_version(psi),
        log_norm = log_norm(psi),
        simple = iteration.observables,
        score = iteration.ctm_score === nothing ? iteration.simple_score : iteration.ctm_score,
    )
    open(path, "w") do io
        JSON3.write(io, payload)
        write(io, '\n')
    end
    return path
end
```

- [ ] **Step 4: Thread store through `scarfinder!`**

Add keyword:

```julia
candidate_store::CandidateStore = NoCandidateStore()
```

After pushing each iteration:

```julia
store_candidate!(candidate_store, psi, iterations[end])
```

- [ ] **Step 5: Run tests**

```bash
julia --project=. test/runtests.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/ScarFinder.jl src/SquarePXPDynamics.jl test/test_scarfinder.jl
git commit -m "feat: persist ScarFinder candidate metadata"
```

## Task 8: Reverse Evolution And Trotter Convergence Checks

**Files:**
- Modify: `src/IPEPSEvolution.jl`
- Modify: `src/PXPValidation.jl`
- Test: `test/test_ipeps_evolution.jl`
- Test: `test/test_pxp_validation.jl`

- [ ] **Step 1: Write reverse API tests**

```julia
@testset "reverse real-time evolution returns close in D1 product limit" begin
    cell = PeriodicSquareUnitCell(10, 10)
    psi = product_square_ipeps(cell; state = :down, maxdim = 1)
    before = measure_simple(psi)
    params = TrotterParams(0.01, 1, :real, 1, 1e-12; schedule = :serial)

    evolve!(psi, 0.01; params)
    reverse_evolve!(psi, 0.01; params)
    after = measure_simple(psi)

    @test after.density ≈ before.density atol = 1e-10
    @test after.blockade_violation ≈ before.blockade_violation atol = 1e-10
end
```

- [ ] **Step 2: Add reverse helper**

In `src/IPEPSEvolution.jl`:

```julia
export reverse_evolve!

"""Evolve backward in real time by applying the same model protocol with negative local steps."""
function reverse_evolve!(
    psi::SquareIPEPSState,
    total_time::Real;
    params::TrotterParams,
    protocol = nothing,
)
    params.evolution === :real ||
        throw(ArgumentError("reverse_evolve! is defined only for real-time evolution"))
    reverse_params = TrotterParams(
        params.dt,
        params.order,
        :real,
        params.maxdim,
        params.cutoff,
        params.split_order;
        schedule = params.schedule,
    )
    return _evolve_with_signed_time!(psi, -Float64(total_time), reverse_params, protocol)
end
```

If `_evolve_with_signed_time!` does not exist, refactor the private evolution loop so the public `evolve!` passes positive signed time and `reverse_evolve!` passes negative signed time, while preserving `TrotterParams.dt > 0`.

- [ ] **Step 3: Add validation summary**

In `src/PXPValidation.jl`, add `PXPReversibilityReport` with before/after simple observables and absolute density/blockade/energy drift.

- [ ] **Step 4: Run tests**

```bash
julia --project=. test/runtests.jl test_ipeps_evolution.jl test_pxp_validation.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/IPEPSEvolution.jl src/PXPValidation.jl src/SquarePXPDynamics.jl test/test_ipeps_evolution.jl test/test_pxp_validation.jl
git commit -m "feat: add PXP reverse-evolution validation"
```

## Task 9: Projection Semantics Clarification

**Files:**
- Modify: `src/SquarePXP.jl`
- Modify: `src/StarModels.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_square_pxp.jl`

- [ ] **Step 1: Write projection semantics tests**

```julia
@testset "projected PXP gate documents left and sandwich projection" begin
    dt = 0.01
    P = square_star_blockade_projector()
    U = square_pxp_gate(dt; evolution = :real)

    @test projected_square_pxp_gate(dt; evolution = :real) ≈ P * U
    @test projected_square_pxp_gate(dt; evolution = :real, projection = :sandwich) ≈ P * U * P
    @test_throws ArgumentError projected_square_pxp_gate(dt; projection = :bad)
end
```

- [ ] **Step 2: Implement keyword**

In `src/SquarePXP.jl`:

```julia
function projected_square_pxp_gate(step::Real; evolution::Symbol = :real, projection::Symbol = :left)
    P = square_star_blockade_projector()
    U = square_pxp_gate(step; evolution)
    projection === :left && return P * U
    projection === :sandwich && return P * U * P
    throw(ArgumentError("projection must be :left or :sandwich"))
end
```

Update docstring to say `:left` preserves the historical assumption that the input is already constrained, while `:sandwich` is explicit constrained-sector action for raw local vectors.

- [ ] **Step 3: Thread model keyword only if needed**

If `PXPStarModel` should expose this, add `projection::Symbol` to the model with default `:left`. If this touches many tests, keep model default unchanged and expose only the dense helper in this task.

- [ ] **Step 4: Run tests**

```bash
julia --project=. test/runtests.jl test_square_pxp.jl test_square_ipeps_s2.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/SquarePXP.jl src/StarModels.jl src/SquarePXPDynamics.jl test/test_square_pxp.jl
git commit -m "docs: clarify PXP projection semantics"
```

## Task 10: CTM-Aware Update Boundary

**Files:**
- Modify: `src/CTMGaugeReadiness.jl`
- Modify: `src/StarSimpleUpdate.jl`
- Modify: `src/SquarePXPDynamics.jl`
- Test: `test/test_ctm_gauge_readiness.jl`

- [ ] **Step 1: Add compatibility wrapper tests**

```julia
@testset "PEPSKit private full-update helpers are wrapped behind compatibility API" begin
    @test SquarePXPDynamics.CTMGaugeReadinessModule.pepskit_private_full_update_available() isa Bool
end
```

- [ ] **Step 2: Wrap private PEPSKit calls**

In `src/CTMGaugeReadiness.jl`, add:

```julia
"""Return true when the PEPSKit private helpers required by CTM gauge conditioning are available."""
function pepskit_private_full_update_available()
    return hasproperty(PEPSKit, :_qr_bond) &&
           hasproperty(PEPSKit, :bondenv_fu) &&
           hasproperty(PEPSKit, :_fixgauge_benvXY)
end
```

Before calling any private helper, add:

```julia
pepskit_private_full_update_available() ||
    throw(ArgumentError("installed PEPSKit does not expose required full-update helper functions"))
```

- [ ] **Step 3: Add CTM-measured simple-update acceptance helper**

In `src/StarSimpleUpdate.jl`, add a small public wrapper only if ScarFinder needs it directly:

```julia
"""Apply simple update, then measure with a trusted CTM backend and restore on trust rejection."""
function project_star_ctm_checked!(psi, center, step; backend, kwargs...)
    before = copy_state(psi)
    info = project_star!(psi, center, step; kwargs...)
    trusted = measure_scarfinder(psi, backend)
    if trusted isa TrustedCTMMeasurement && !trusted.trust.trusted
        _replace_state!(psi, before)
    end
    return info, trusted
end
```

If this introduces circular imports, do not add it here; keep CTM-measured acceptance in ScarFinder and document full update as future infrastructure.

- [ ] **Step 4: Run focused tests**

```bash
julia --project=. test/runtests.jl test_ctm_gauge_readiness.jl test_scarfinder.jl test_public_docs.jl
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/CTMGaugeReadiness.jl src/StarSimpleUpdate.jl src/SquarePXPDynamics.jl test/test_ctm_gauge_readiness.jl
git commit -m "refactor: wrap PEPSKit CTM update compatibility"
```

## Task 11: README And Final Verification

**Files:**
- Modify: `README.md`
- Modify: `memory/mid_term/decision_log.md`

- [ ] **Step 1: Update README status**

Update the ScarFinder section so it states:

```markdown
ScarFinder now supports two measurement modes: `SimpleBackend()` for fast
development diagnostics and `TrustedCTMBackend(...)` for finite-chi CTM-backed
candidate ranking. Physics-facing runs should use `TrustedCTMBackend`,
`CompositeObjective`, and convergence reports from `validate_pxp_convergence`.
Simple/local diagnostics remain useful for smoke tests and regression checks,
but not for final claims.
```

- [ ] **Step 2: Record architecture decision**

Append to `memory/mid_term/decision_log.md`:

```markdown
## 2026-05-16 - ScarFinder Uses Trusted Measurement Backends

Decision:

Promote ScarFinder from callback-based optional CTM diagnostics to explicit
measurement backends and physics objectives. `SimpleBackend` remains the fast
development path, while `TrustedCTMBackend` is the production ranking path.

Reason:

The GPT review correctly identified that simple/local diagnostics cannot be
the default evidence layer for ScarFinder physics claims.

Consequences:

Candidate ranking can now require finite-chi CTM trust, store trust policy and
drift metadata, and be audited alongside ED and convergence reports.

Source:

`src/ScarFinder.jl`; `src/PXPValidation.jl`; this implementation plan.

Status: active
```

- [ ] **Step 3: Run full verification**

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
git diff --check
```

Expected: full package tests PASS and `git diff --check` prints no output.

- [ ] **Step 4: Request code review**

Use `superpowers:requesting-code-review` and ask the reviewer to check:

- trusted CTM ranking cannot accidentally rank untrusted CTM candidates when `require_ctm_trusted = true`
- objective score signs match the intended convention: smaller is better, stronger revival lowers score
- deterministic CTMRG seeding does not break unseeded PEPSKit behavior
- JSON candidate store does not claim to serialize full tensor states unless full state persistence is added
- projection semantics preserve the current default `P * U`

- [ ] **Step 5: Commit docs**

```bash
git add README.md memory/mid_term/decision_log.md
git commit -m "docs: document trusted ScarFinder workflow"
```

## Final Definition Of Done

- [ ] `scarfinder!` accepts `measurement = TrustedCTMBackend(...)`, `objective = CompositeObjective(...)`, and `require_trusted_ctm = true`.
- [ ] `rank_scarfinder_candidates(...; diagnostics = :ctm, require_ctm_trusted = true)` excludes untrusted finite-chi CTM records.
- [ ] CTM summaries include sublattice imbalance and checkerboard structure-factor diagnostics.
- [ ] CTMRG parameters can record a seed and serialized artifacts include it.
- [ ] `validate_pxp_convergence` produces a machine-readable dt/D/cutoff/chi report.
- [ ] Candidate metadata is persisted per iteration.
- [ ] Reverse real-time validation exists for short D=1/D=2 smoke checks.
- [ ] `projected_square_pxp_gate` documents current `P * U` behavior and exposes `P * U * P`.
- [ ] PEPSKit private helper usage is isolated behind compatibility checks.
- [ ] README states the new trusted ScarFinder workflow and remaining limits.
- [ ] `julia --project=. -e 'using Pkg; Pkg.test()'` passes.
- [ ] `git diff --check` passes.

## Self-Review

- Spec coverage: The plan covers all remaining GPT comments except a full ALS/full-update solver as a production algorithm. That item is deliberately scoped as CTM-measured simple update plus PEPSKit compatibility boundaries first, matching the review's recommended option 1 before option 2.
- Placeholder scan: No task uses unspecified implementation language; every task includes concrete tests, code shapes, commands, and expected outcomes.
- Type consistency: New public names are consistently `MeasurementBackend`, `SimpleBackend`, `TrustedCTMBackend`, `ScarFinderObjective`, `CompositeObjective`, `PXPConvergenceConfig`, and `PXPConvergenceReport`.
