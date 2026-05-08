module ScarFinder

using ..States: AbstractUnitCell, TriangularIPEPS, OneSiteUnitCell, product_ipeps, random_ipeps,
                unit_cell_representatives
using ..Evolution: ProjectedPXPStepDiagnostics, run_projected_pxp!
using ..Observables: mean_blockade_violation

export ScarFinderConfig, ScarCandidate, scar_search, rank_candidates

struct ScarFinderConfig
    dt::Float64
    projection_interval::Int
    niterations::Int
    maxdim::Int
    cutoff::Float64
    unitcell::AbstractUnitCell
    seed_count::Int
    blockade_tolerance::Float64

    function ScarFinderConfig(dt::Real,
                              projection_interval::Integer,
                              niterations::Integer,
                              maxdim::Integer,
                              cutoff::Real,
                              unitcell::AbstractUnitCell,
                              seed_count::Integer,
                              blockade_tolerance::Real)
        dt > 0 || throw(ArgumentError("dt must be positive"))
        projection_interval >= 1 || throw(ArgumentError("projection_interval must be >= 1"))
        niterations >= 0 || throw(ArgumentError("niterations must be nonnegative"))
        maxdim >= 1 || throw(ArgumentError("maxdim must be >= 1"))
        cutoff >= 0 || throw(ArgumentError("cutoff must be nonnegative"))
        seed_count >= 1 || throw(ArgumentError("seed_count must be >= 1"))
        blockade_tolerance >= 0 || throw(ArgumentError("blockade_tolerance must be nonnegative"))
        return new(Float64(dt), Int(projection_interval), Int(niterations),
                   Int(maxdim), Float64(cutoff), unitcell, Int(seed_count),
                   Float64(blockade_tolerance))
    end
end

struct ScarCandidate
    seed_index::Int
    seed_kind::Symbol
    state::TriangularIPEPS
    diagnostics::Vector{ProjectedPXPStepDiagnostics}
    discarded_weight::Float64
    blockade_violation::Float64
    entanglement_proxy::Float64
    score::Float64
    accepted::Bool
end

function scar_search(config::ScarFinderConfig; seed::Integer = 0)
    config.maxdim == 1 || throw(ArgumentError(
        "scar_search currently supports only D=1 product-state searches; " *
        "fixed-D PEPS ScarFinder requires a general D>1 Simple Update backend"
    ))
    candidates = ScarCandidate[]
    for seed_index in 1:config.seed_count
        kind, state = _seed_state(config, seed_index, seed)
        diagnostics = ProjectedPXPStepDiagnostics[]
        failed_projection = false
        for _ in 1:config.niterations
            try
                append!(diagnostics, run_projected_pxp!(
                    state, config.dt, config.projection_interval;
                    order = :second, maxdim = config.maxdim, cutoff = config.cutoff,
                    evolution = :real,
                ))
            catch err
                err isa ArgumentError || rethrow()
                failed_projection = true
                break
            end
        end
        discarded = sum((sum(d.discarded_weights) for d in diagnostics); init = 0.0)
        blockade = mean_blockade_violation(state, unit_cell_representatives(config.unitcell))
        entropy = _lambda_entropy_proxy(state)
        score = discarded + blockade + entropy + (failed_projection ? 1.0e6 : 0.0)
        accepted = !failed_projection && blockade <= config.blockade_tolerance
        push!(candidates, ScarCandidate(seed_index, kind, state, diagnostics,
                                        discarded, blockade, entropy, score,
                                        accepted))
    end
    return rank_candidates(candidates)
end

function rank_candidates(candidates)
    return sort(collect(candidates); by = c -> (c.score, c.blockade_violation,
                                                c.entanglement_proxy, c.seed_index))
end

function _seed_state(config::ScarFinderConfig, seed_index::Int, seed::Integer)
    if seed_index == 1
        return :product_down, product_ipeps(config.unitcell, :down; D = config.maxdim)
    elseif seed_index == 2
        return :product_up, product_ipeps(config.unitcell, :up; D = config.maxdim)
    else
        return :random, random_ipeps(config.unitcell, config.maxdim; seed = seed + seed_index)
    end
end

function _lambda_entropy_proxy(state::TriangularIPEPS)
    entropies = Float64[]
    seen = Set{UInt}()
    for λ in values(state.lambdas)
        objectid(λ) in seen && continue
        push!(seen, objectid(λ))
        weights = abs2.(λ)
        total = sum(weights)
        if total == 0
            push!(entropies, 0.0)
        else
            p = weights ./ total
            push!(entropies, -sum(pi > 0 ? pi * log(pi) : 0.0 for pi in p))
        end
    end
    return isempty(entropies) ? 0.0 : sum(entropies) / length(entropies)
end

end
