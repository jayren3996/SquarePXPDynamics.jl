module ScarFinderAudit

using JSON3
using Statistics: mean, std

using ..Internals: _csv_value
using ..SquareUnitCells: PeriodicSquareUnitCell
using ..SquareIPEPS: SquareIPEPSState, product_square_ipeps, checkerboard_square_ipeps
using ..IPEPSEvolution: TrotterParams
using ..PEPSKitMeasurements: PEPSKitCTMRGParams
using ..CTMTrust: CTMTrustPolicy
using ..PXPValidation: validate_pxp_reversibility, PXPReversibilityReport
using ..ScarFinder:
    ScarFinderParams,
    ScarFinderResult,
    ScarFinderIteration,
    ScarFinderCandidateScore,
    SimpleBackend,
    TrustedCTMBackend,
    MeasurementBackend,
    CompositeObjective,
    ScarFinderObjective,
    CandidateStore,
    NoCandidateStore,
    scarfinder!

export ScarFinderAuditConfig,
    ScarFinderAuditRow,
    ScarFinderAuditStability,
    ScarFinderAuditReport,
    run_scarfinder_audit,
    write_scarfinder_audit_json,
    write_scarfinder_audit_csv

"""
    ScarFinderAuditConfig(; kwargs...)

Grid specification for [`run_scarfinder_audit`](@ref). Each combination of
`(dt, target_maxdim, cutoff)` is run once with the simple-update backend.
Each combination of `(dt, target_maxdim, cutoff, chi)` is additionally run
with a trusted-CTM backend when `chi_values` is nonempty.

`evolve_maxdim` is the **fixed** simple-update bond dimension during
`evolve!` (the physical-evolution accuracy axis — held constant across the
audit). Each `target_maxdim` in `target_maxdim_values` is the
scarfinder-level compression target after each iteration's `evolve!`
(the scar-candidate axis the audit sweeps). All `target_maxdim` values must
satisfy `target_maxdim <= evolve_maxdim`.

`projection_time` is the **per-iteration** evolution duration applied by
`scarfinder!`; total evolution over a row is `projection_time * iterations`.
States are constructed fresh from `initial_state` at the row's
`target_maxdim`, then naturally grow up to `evolve_maxdim` during `evolve!`.
"""
struct ScarFinderAuditConfig
    cell_Lx::Int
    cell_Ly::Int
    initial_state::Symbol
    projection_time::Float64
    iterations::Int
    dt_values::Vector{Float64}
    evolve_maxdim::Int
    target_maxdim_values::Vector{Int}
    cutoff_values::Vector{Float64}
    chi_values::Vector{Int}
    schedule::Symbol
    order::Int
    blockade_tol::Float64
    truncation_tol::Float64
    entropy_cap::Float64
    compression_cutoff::Float64
    ctm_tol::Float64
    ctm_maxiter::Int
    ctm_seed::Union{Nothing,Int}
    ctm_trust_policy::CTMTrustPolicy
    require_trusted_ctm::Bool
    top_k::Int

    function ScarFinderAuditConfig(;
        cell_Lx::Integer,
        cell_Ly::Integer,
        initial_state::Symbol = :down,
        projection_time::Real,
        iterations::Integer,
        dt_values,
        evolve_maxdim::Integer,
        target_maxdim_values,
        cutoff_values = [1e-12],
        chi_values = Int[],
        schedule::Symbol = :serial,
        order::Integer = 1,
        blockade_tol::Real = Inf,
        truncation_tol::Real = Inf,
        entropy_cap::Real = Inf,
        compression_cutoff::Real = 1e-12,
        ctm_tol::Real = 1e-8,
        ctm_maxiter::Integer = 100,
        ctm_seed::Union{Nothing,Integer} = 0,
        ctm_trust_policy::CTMTrustPolicy = CTMTrustPolicy(),
        require_trusted_ctm::Bool = false,
        top_k::Integer = 3,
    )
        cell_Lx >= 2 || throw(ArgumentError("cell_Lx must be at least 2"))
        cell_Ly >= 2 || throw(ArgumentError("cell_Ly must be at least 2"))
        iterations >= 1 || throw(ArgumentError("iterations must be at least 1"))
        projection_time >= 0 || throw(ArgumentError("projection_time must be nonnegative"))
        order in (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        schedule in (:serial, :five_color) ||
            throw(ArgumentError("schedule must be :serial or :five_color"))
        top_k >= 1 || throw(ArgumentError("top_k must be at least 1"))
        evolve_maxdim >= 1 || throw(ArgumentError("evolve_maxdim must be at least 1"))
        compression_cutoff >= 0 && isfinite(compression_cutoff) ||
            throw(ArgumentError("compression_cutoff must be finite and nonnegative"))

        dt_grid = [Float64(dt) for dt in dt_values]
        target_grid = [Int(D) for D in target_maxdim_values]
        cutoff_grid = [Float64(cutoff) for cutoff in cutoff_values]
        chi_grid = [Int(chi) for chi in chi_values]

        isempty(dt_grid) && throw(ArgumentError("dt_values must be nonempty"))
        isempty(target_grid) && throw(ArgumentError("target_maxdim_values must be nonempty"))
        isempty(cutoff_grid) && throw(ArgumentError("cutoff_values must be nonempty"))

        all(dt -> dt > 0 && isfinite(dt), dt_grid) ||
            throw(ArgumentError("dt_values must be positive and finite"))
        all(D -> D >= 1, target_grid) ||
            throw(ArgumentError("target_maxdim_values must be at least 1"))
        all(D -> D <= evolve_maxdim, target_grid) ||
            throw(ArgumentError("target_maxdim_values must be <= evolve_maxdim"))
        all(c -> c >= 0 && isfinite(c), cutoff_grid) ||
            throw(ArgumentError("cutoff_values must be nonnegative and finite"))
        all(chi -> chi >= 1, chi_grid) ||
            throw(ArgumentError("chi_values must be at least 1"))

        require_trusted_ctm && isempty(chi_grid) &&
            throw(ArgumentError("require_trusted_ctm requires nonempty chi_values"))

        return new(
            Int(cell_Lx),
            Int(cell_Ly),
            initial_state,
            Float64(projection_time),
            Int(iterations),
            dt_grid,
            Int(evolve_maxdim),
            target_grid,
            cutoff_grid,
            chi_grid,
            schedule,
            Int(order),
            Float64(blockade_tol),
            Float64(truncation_tol),
            Float64(entropy_cap),
            Float64(compression_cutoff),
            Float64(ctm_tol),
            Int(ctm_maxiter),
            ctm_seed === nothing ? nothing : Int(ctm_seed),
            ctm_trust_policy,
            require_trusted_ctm,
            Int(top_k),
        )
    end
end

"""
    ScarFinderAuditRow

One row of a [`ScarFinderAuditReport`](@ref). `backend` is `:simple` for
simple-update rows and `:ctm` for trusted-CTM rows. `per_iteration_scores`
contains the per-iteration `score` field for the chosen backend, with `NaN`
in slots where the iteration was rejected and produced no score. `chi` is
`nothing` for simple-backend rows. `target_maxdim` is the scarfinder
compression cap for this row; `evolve_maxdim` is the simulation accuracy
shared across the audit.

`max_truncerr_used` is the largest `evolution.max_truncerr` observed across
this row's scarfinder! iterations (acceptance criterion 6 — iPEPS truncation
budget). `max_compression_truncerr` is the largest per-bond compression
truncation error across iterations (the scar-quality signal). `reversibility_*_drift`
fields hold the absolute forward/reverse drift in density, blockade violation,
and PXP energy density for the `(dt, target_maxdim, cutoff)` combination,
computed once and shared across this row's chi values (acceptance criterion 5).
All extra fields are `nothing` when the row errored before they could be populated.
"""
struct ScarFinderAuditRow
    dt::Float64
    target_maxdim::Int
    evolve_maxdim::Int
    cutoff::Float64
    chi::Union{Nothing,Int}
    backend::Symbol
    accepted::Int
    rejected::Int
    iterations_run::Int
    best_iteration::Union{Nothing,Int}
    best_score::Union{Nothing,Float64}
    top_k_iteration_indices::Vector{Int}
    per_iteration_scores::Vector{Float64}
    ctm_trusted_count::Union{Nothing,Int}
    max_truncerr_used::Union{Nothing,Float64}
    max_compression_truncerr::Union{Nothing,Float64}
    mean_compression_truncerr::Union{Nothing,Float64}
    max_compression_density_shift::Union{Nothing,Float64}
    reversibility_density_drift::Union{Nothing,Float64}
    reversibility_blockade_drift::Union{Nothing,Float64}
    reversibility_energy_drift::Union{Nothing,Float64}
    elapsed_seconds::Float64
    ok::Bool
    error::Union{Nothing,String}
end

"""
    ScarFinderAuditStability

Cross-row summary diagnostics from [`run_scarfinder_audit`](@ref). Each metric
is computed three ways: across all successful rows (`*_all`), across only
simple-backend rows (`*_simple`), and across only CTM-backend rows (`*_ctm`).
Mixing backends compares incomparable score conventions; the per-backend
breakdown is the authoritative measure for ranking stability.

`best_iteration_agreement` is the fraction of rows whose best-iteration index
equals the modal best-iteration index in the relevant subset.
`score_range` and `score_cv` summarize `best_score` variation across rows in
the subset. `top_k_jaccard_min` is the smallest pairwise Jaccard overlap
between rows' `top_k_iteration_indices`. All fields are `nothing` when fewer
than two rows are available in the relevant subset.
"""
struct ScarFinderAuditStability
    n_rows::Int
    n_successful::Int
    n_simple::Int
    n_ctm::Int
    best_iteration_agreement_all::Union{Nothing,Float64}
    score_range_all::Union{Nothing,Float64}
    score_cv_all::Union{Nothing,Float64}
    top_k_jaccard_min_all::Union{Nothing,Float64}
    best_iteration_agreement_simple::Union{Nothing,Float64}
    score_range_simple::Union{Nothing,Float64}
    score_cv_simple::Union{Nothing,Float64}
    top_k_jaccard_min_simple::Union{Nothing,Float64}
    best_iteration_agreement_ctm::Union{Nothing,Float64}
    score_range_ctm::Union{Nothing,Float64}
    score_cv_ctm::Union{Nothing,Float64}
    top_k_jaccard_min_ctm::Union{Nothing,Float64}
end

"""
    ScarFinderAuditReport(config, rows, stability)

Full result of [`run_scarfinder_audit`](@ref).
"""
struct ScarFinderAuditReport
    config::ScarFinderAuditConfig
    rows::Vector{ScarFinderAuditRow}
    stability::ScarFinderAuditStability
end

function _build_state(cfg::ScarFinderAuditConfig, target_maxdim::Int)
    cell = PeriodicSquareUnitCell(cfg.cell_Lx, cfg.cell_Ly)
    if cfg.initial_state === :checkerboard_even
        return checkerboard_square_ipeps(cell; excited_on = :even, maxdim = target_maxdim)
    elseif cfg.initial_state === :checkerboard_odd
        return checkerboard_square_ipeps(cell; excited_on = :odd, maxdim = target_maxdim)
    else
        return product_square_ipeps(cell; state = cfg.initial_state, maxdim = target_maxdim)
    end
end

function _scarfinder_params(cfg::ScarFinderAuditConfig, dt::Float64, target_maxdim::Int,
                            cutoff::Float64)
    trotter = TrotterParams(
        dt, cfg.order, :real, cfg.evolve_maxdim, cutoff; schedule = cfg.schedule,
    )
    return ScarFinderParams(
        cfg.projection_time,
        trotter,
        cfg.iterations,
        cfg.blockade_tol,
        cfg.truncation_tol,
        cfg.entropy_cap,
        false;
        target_maxdim = target_maxdim,
        compression_cutoff = cfg.compression_cutoff,
    )
end

function _score_for(iteration::ScarFinderIteration, backend::Symbol)
    if backend === :simple
        return iteration.simple_score
    elseif backend === :ctm
        return iteration.ctm_score
    else
        throw(ArgumentError("unknown audit backend $backend"))
    end
end

function _per_iteration_scores(result::ScarFinderResult, backend::Symbol, planned::Int)
    scores = fill(NaN, planned)
    for iteration in result.iterations
        score = _score_for(iteration, backend)
        score === nothing && continue
        scores[iteration.iteration] = score.score
    end
    return scores
end

function _top_k_iterations(scores::Vector{Float64}, k::Int)
    indices = collect(eachindex(scores))
    finite = filter(i -> isfinite(scores[i]), indices)
    isempty(finite) && return Int[]
    sort!(finite; by = i -> scores[i], rev = true)
    return finite[1:min(k, length(finite))]
end

function _best_finite(scores::Vector{Float64})
    finite = [i for i in eachindex(scores) if isfinite(scores[i])]
    isempty(finite) && return (nothing, nothing)
    idx = argmax(scores[finite])
    chosen = finite[idx]
    return (chosen, scores[chosen])
end

function _count_trusted_ctm(result::ScarFinderResult)
    return count(iteration -> iteration.ctm_score !== nothing, result.iterations)
end

function _max_truncerr(result::ScarFinderResult)
    isempty(result.iterations) && return nothing
    raw = [iteration.evolution.max_truncerr for iteration in result.iterations]
    finite = filter(isfinite, raw)
    isempty(finite) && return nothing
    return maximum(finite)
end

function _reversibility_report(cfg::ScarFinderAuditConfig, dt::Float64,
                               target_maxdim::Int, cutoff::Float64)
    cell = PeriodicSquareUnitCell(cfg.cell_Lx, cfg.cell_Ly)
    psi = product_square_ipeps(cell; state = cfg.initial_state, maxdim = target_maxdim)
    trotter = TrotterParams(
        dt, cfg.order, :real, cfg.evolve_maxdim, cutoff; schedule = cfg.schedule,
    )
    total = cfg.projection_time * cfg.iterations
    return validate_pxp_reversibility(psi, total; params = trotter)
end

function _run_simple_row(cfg::ScarFinderAuditConfig, dt::Float64, target_maxdim::Int,
                         cutoff::Float64, objective::ScarFinderObjective,
                         rev::Union{Nothing,PXPReversibilityReport})
    psi = _build_state(cfg, target_maxdim)
    params = _scarfinder_params(cfg, dt, target_maxdim, cutoff)
    elapsed = NaN
    err = nothing
    result = nothing
    try
        elapsed = @elapsed begin
            result = scarfinder!(psi, params; objective, measurement = SimpleBackend())
        end
    catch e
        err = sprint(showerror, e)
    end
    return _row_from_result(
        result, :simple, dt, target_maxdim, cutoff, nothing, cfg, elapsed, err, rev,
    )
end

function _trust_chi_window(target_chi::Int, min_points::Int)
    points = max(min_points, 2)
    raw = [target_chi - 4 * i for i in (points - 1):-1:0]
    cleaned = max.(raw, 2)
    return unique(sort(cleaned))
end

function _run_ctm_row(cfg::ScarFinderAuditConfig, dt::Float64, target_maxdim::Int,
                      cutoff::Float64, chi::Int, objective::ScarFinderObjective,
                      rev::Union{Nothing,PXPReversibilityReport})
    psi = _build_state(cfg, target_maxdim)
    params = _scarfinder_params(cfg, dt, target_maxdim, cutoff)
    window_chis = _trust_chi_window(chi, cfg.ctm_trust_policy.min_points)
    ctm_window = Tuple(
        PEPSKitCTMRGParams(c, cfg.ctm_tol, cfg.ctm_maxiter, 0; seed = cfg.ctm_seed)
        for c in window_chis
    )
    backend = TrustedCTMBackend(ctm_window, cfg.ctm_trust_policy)
    elapsed = NaN
    err = nothing
    result = nothing
    try
        elapsed = @elapsed begin
            result = scarfinder!(
                psi,
                params;
                objective,
                measurement = backend,
                ctm_every = 1,
                require_trusted_ctm = cfg.require_trusted_ctm,
            )
        end
    catch e
        err = sprint(showerror, e)
    end
    return _row_from_result(
        result, :ctm, dt, target_maxdim, cutoff, chi, cfg, elapsed, err, rev,
    )
end

function _compression_stats(result::ScarFinderResult)
    compressions = [it.compression for it in result.iterations if it.compression !== nothing]
    isempty(compressions) && return (nothing, nothing, nothing)
    max_trunc = maximum(c.info.max_truncerr for c in compressions)
    mean_trunc = sum(c.info.mean_truncerr for c in compressions) / length(compressions)
    max_density_shift = maximum(c.density_shift for c in compressions)
    return (max_trunc, mean_trunc, max_density_shift)
end

function _row_from_result(result::Union{Nothing,ScarFinderResult}, backend::Symbol,
                          dt::Float64, target_maxdim::Int, cutoff::Float64,
                          chi::Union{Nothing,Int}, cfg::ScarFinderAuditConfig,
                          elapsed::Float64, err::Union{Nothing,String},
                          rev::Union{Nothing,PXPReversibilityReport})
    rev_density = rev === nothing ? nothing : rev.density_drift
    rev_blockade = rev === nothing ? nothing : rev.blockade_drift
    rev_energy = rev === nothing ? nothing : rev.energy_drift
    if result === nothing
        return ScarFinderAuditRow(
            dt, target_maxdim, cfg.evolve_maxdim, cutoff, chi, backend, 0, 0, 0,
            nothing, nothing, Int[], Float64[],
            backend === :ctm ? 0 : nothing,
            nothing, nothing, nothing, nothing,
            rev_density, rev_blockade, rev_energy,
            elapsed, false, err,
        )
    end
    scores = _per_iteration_scores(result, backend, cfg.iterations)
    best_idx, best_score = _best_finite(scores)
    top_k = _top_k_iterations(scores, cfg.top_k)
    max_comp_trunc, mean_comp_trunc, max_comp_density_shift = _compression_stats(result)
    return ScarFinderAuditRow(
        dt, target_maxdim, cfg.evolve_maxdim, cutoff, chi, backend,
        result.accepted_iterations,
        result.rejected_iterations,
        length(result.iterations),
        best_idx,
        best_score,
        top_k,
        scores,
        backend === :ctm ? _count_trusted_ctm(result) : nothing,
        _max_truncerr(result),
        max_comp_trunc,
        mean_comp_trunc,
        max_comp_density_shift,
        rev_density,
        rev_blockade,
        rev_energy,
        elapsed,
        true,
        err,
    )
end

function _subset_metrics(rows::Vector{ScarFinderAuditRow})
    if length(rows) < 2
        return (nothing, nothing, nothing, nothing)
    end
    best_iters = [row.best_iteration::Int for row in rows]
    modal = _modal_value(best_iters)
    agreement = count(==(modal), best_iters) / length(best_iters)
    best_scores = [row.best_score::Float64 for row in rows]
    score_range = maximum(best_scores) - minimum(best_scores)
    avg = mean(best_scores)
    score_cv = abs(avg) > eps() ? std(best_scores) / abs(avg) : nothing
    jaccard_min =
        _pairwise_jaccard_min([row.top_k_iteration_indices for row in rows])
    return (agreement, score_range, score_cv, jaccard_min)
end

function _stability(rows::Vector{ScarFinderAuditRow})
    successful = filter(row -> row.ok && row.best_score !== nothing, rows)
    simple_rows = filter(row -> row.backend === :simple, successful)
    ctm_rows = filter(row -> row.backend === :ctm, successful)

    all_metrics = _subset_metrics(successful)
    simple_metrics = _subset_metrics(simple_rows)
    ctm_metrics = _subset_metrics(ctm_rows)

    return ScarFinderAuditStability(
        length(rows),
        length(successful),
        length(simple_rows),
        length(ctm_rows),
        all_metrics...,
        simple_metrics...,
        ctm_metrics...,
    )
end

function _modal_value(values)
    counts = Dict{eltype(values),Int}()
    for v in values
        counts[v] = get(counts, v, 0) + 1
    end
    return argmax(counts)
end

function _pairwise_jaccard_min(top_k_lists::Vector{Vector{Int}})
    isempty(top_k_lists) && return nothing
    length(top_k_lists) < 2 && return nothing
    min_overlap = 1.0
    for i in 1:(length(top_k_lists) - 1), j in (i + 1):length(top_k_lists)
        a = Set(top_k_lists[i])
        b = Set(top_k_lists[j])
        union_size = length(union(a, b))
        union_size == 0 && continue
        overlap = length(intersect(a, b)) / union_size
        min_overlap = min(min_overlap, overlap)
    end
    return min_overlap
end

"""
    run_scarfinder_audit(config; objective = CompositeObjective())

Sweep `(dt, target_maxdim, cutoff)` simple-backend ScarFinder runs and, when
`config.chi_values` is nonempty, `(dt, target_maxdim, cutoff, chi)` trusted-CTM
runs. Each run evolves at `config.evolve_maxdim` and compresses to
`target_maxdim` after every iteration. Returns a [`ScarFinderAuditReport`](@ref)
with per-row diagnostics and a cross-row stability summary.
"""
function run_scarfinder_audit(
    config::ScarFinderAuditConfig;
    objective::ScarFinderObjective = CompositeObjective(),
)
    rows = ScarFinderAuditRow[]
    for dt in config.dt_values,
        target_maxdim in config.target_maxdim_values,
        cutoff in config.cutoff_values

        rev = try
            _reversibility_report(config, dt, target_maxdim, cutoff)
        catch err
            @warn "reversibility check failed" dt target_maxdim cutoff exception = err
            nothing
        end
        push!(rows, _run_simple_row(config, dt, target_maxdim, cutoff, objective, rev))
        for chi in config.chi_values
            push!(
                rows,
                _run_ctm_row(config, dt, target_maxdim, cutoff, chi, objective, rev),
            )
        end
    end
    return ScarFinderAuditReport(config, rows, _stability(rows))
end

const _AUDIT_CSV_COLUMNS = (
    :dt, :target_maxdim, :evolve_maxdim, :cutoff, :chi, :backend,
    :accepted, :rejected, :iterations_run,
    :best_iteration, :best_score,
    :top_k_iteration_indices, :ctm_trusted_count,
    :max_truncerr_used,
    :max_compression_truncerr,
    :mean_compression_truncerr,
    :max_compression_density_shift,
    :reversibility_density_drift,
    :reversibility_blockade_drift,
    :reversibility_energy_drift,
    :elapsed_seconds, :ok, :error,
)

function _csv_format(value)
    return value === nothing ? "" : _csv_value(value)
end

function _format_top_k(indices::Vector{Int})
    isempty(indices) && return ""
    return "\"" * join(string.(indices), "|") * "\""
end

"""
    write_scarfinder_audit_csv(report, path)

Write per-row audit results to `path` as CSV. The `top_k_iteration_indices`
column is pipe-delimited inside quotes; `per_iteration_scores` is omitted from
CSV — use JSON output to read it.
"""
function write_scarfinder_audit_csv(report::ScarFinderAuditReport, path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(string.(_AUDIT_CSV_COLUMNS), ","))
        for row in report.rows
            fields = String[]
            for col in _AUDIT_CSV_COLUMNS
                value = getproperty(row, col)
                push!(fields, col === :top_k_iteration_indices ? _format_top_k(value) : _csv_format(value))
            end
            println(io, join(fields, ","))
        end
    end
    return path
end

_finite_or_nothing(value::Real) = isfinite(value) ? Float64(value) : nothing
_finite_or_nothing(::Nothing) = nothing
_finite_or_nan(value::Real) = isfinite(value) ? Float64(value) : NaN
_finite_or_nan(::Nothing) = nothing

function _config_payload(cfg::ScarFinderAuditConfig)
    return (;
        cell_Lx = cfg.cell_Lx,
        cell_Ly = cfg.cell_Ly,
        initial_state = String(cfg.initial_state),
        projection_time = cfg.projection_time,
        total_evolution_time = cfg.projection_time * cfg.iterations,
        iterations = cfg.iterations,
        dt_values = cfg.dt_values,
        evolve_maxdim = cfg.evolve_maxdim,
        target_maxdim_values = cfg.target_maxdim_values,
        cutoff_values = cfg.cutoff_values,
        chi_values = cfg.chi_values,
        schedule = String(cfg.schedule),
        order = cfg.order,
        blockade_tol = _finite_or_nothing(cfg.blockade_tol),
        truncation_tol = _finite_or_nothing(cfg.truncation_tol),
        entropy_cap = _finite_or_nothing(cfg.entropy_cap),
        compression_cutoff = cfg.compression_cutoff,
        ctm_tol = cfg.ctm_tol,
        ctm_maxiter = cfg.ctm_maxiter,
        ctm_seed = cfg.ctm_seed,
        require_trusted_ctm = cfg.require_trusted_ctm,
        top_k = cfg.top_k,
        ctm_trust_policy = (;
            min_points = cfg.ctm_trust_policy.min_points,
            require_accepted_diagnostics = cfg.ctm_trust_policy.require_accepted_diagnostics,
            max_density_delta = cfg.ctm_trust_policy.max_density_delta,
            max_blockade_delta = cfg.ctm_trust_policy.max_blockade_delta,
            max_energy_delta = cfg.ctm_trust_policy.max_energy_delta,
            max_residual = cfg.ctm_trust_policy.max_residual,
        ),
    )
end

_json_score(value::Float64) = isfinite(value) ? value : nothing
_json_score(::Nothing) = nothing

function _row_payload(row::ScarFinderAuditRow)
    return (;
        dt = row.dt,
        target_maxdim = row.target_maxdim,
        evolve_maxdim = row.evolve_maxdim,
        cutoff = row.cutoff,
        chi = row.chi,
        backend = String(row.backend),
        accepted = row.accepted,
        rejected = row.rejected,
        iterations_run = row.iterations_run,
        best_iteration = row.best_iteration,
        best_score = _json_score(row.best_score),
        top_k_iteration_indices = row.top_k_iteration_indices,
        per_iteration_scores = [_json_score(s) for s in row.per_iteration_scores],
        ctm_trusted_count = row.ctm_trusted_count,
        max_truncerr_used = row.max_truncerr_used,
        max_compression_truncerr = row.max_compression_truncerr,
        mean_compression_truncerr = row.mean_compression_truncerr,
        max_compression_density_shift = row.max_compression_density_shift,
        reversibility_density_drift = row.reversibility_density_drift,
        reversibility_blockade_drift = row.reversibility_blockade_drift,
        reversibility_energy_drift = row.reversibility_energy_drift,
        elapsed_seconds = isfinite(row.elapsed_seconds) ? row.elapsed_seconds : nothing,
        ok = row.ok,
        error = row.error,
    )
end

"""
    write_scarfinder_audit_json(report, path)

Write a full audit report to `path` as JSON, including configuration,
per-row data, and cross-row stability summary.
"""
function write_scarfinder_audit_json(report::ScarFinderAuditReport, path::AbstractString)
    mkpath(dirname(path))
    payload = (;
        config = _config_payload(report.config),
        rows = [_row_payload(row) for row in report.rows],
        stability = (;
            n_rows = report.stability.n_rows,
            n_successful = report.stability.n_successful,
            n_simple = report.stability.n_simple,
            n_ctm = report.stability.n_ctm,
            best_iteration_agreement_all = report.stability.best_iteration_agreement_all,
            score_range_all = report.stability.score_range_all,
            score_cv_all = report.stability.score_cv_all,
            top_k_jaccard_min_all = report.stability.top_k_jaccard_min_all,
            best_iteration_agreement_simple = report.stability.best_iteration_agreement_simple,
            score_range_simple = report.stability.score_range_simple,
            score_cv_simple = report.stability.score_cv_simple,
            top_k_jaccard_min_simple = report.stability.top_k_jaccard_min_simple,
            best_iteration_agreement_ctm = report.stability.best_iteration_agreement_ctm,
            score_range_ctm = report.stability.score_range_ctm,
            score_cv_ctm = report.stability.score_cv_ctm,
            top_k_jaccard_min_ctm = report.stability.top_k_jaccard_min_ctm,
        ),
    )
    open(path, "w") do io
        JSON3.pretty(io, payload)
        println(io)
    end
    return path
end

end # module
