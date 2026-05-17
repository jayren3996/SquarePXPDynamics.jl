module FinitePXPEEDBenchmark

using JSON3
using LinearAlgebra
using SparseArrays

import EDKit

export PXPSquareSpaceGroupBasis
export PXPEEDBenchmarkConfig, PXPEEDSample, PXPEEDBenchmarkResult
export pxp_ed_space_group_basis, pxp_ed_constrained_count, pxp_ed_group_order
export pxp_ed_boundary_condition, pxp_ed_symmetry_sector, pxp_ed_observable_scope
export pxp_ed_reference_label, pxp_ed_site_density_operator, pxp_ed_region_density_operator
export pxp_ed_initial_state, pxp_ed_hamiltonian_operator, sparse_pxp_ed_hamiltonian
export run_pxp_ed_benchmark, write_pxp_ed_benchmark_json

const _MAX_UINT64_SITES = 62

"""
    PXPSquareSpaceGroupBasis

EDKit-compatible basis for periodic square-lattice PXP states reduced by
translations and, by default, the square point group in the fully symmetric
sector.

The physical digit convention follows the rest of the package and EDKit:
digit `0` is the Rydberg/excited `:up` state and digit `1` is `:down`.
"""
struct PXPSquareSpaceGroupBasis <: EDKit.AbstractPermuteBasis
    dgt::Vector{Int}
    I::Vector{Int}
    R::Vector{Float64}
    B::Int
    n::Int
    point_group::Bool
    constrained_count::Int
    perms::Vector{Vector{Int}}
    perm_tables::Vector{Matrix{UInt64}}
    lookup::Dict{Int,Int}
    mask::UInt64
    row_mask::UInt64
end

Base.copy(b::PXPSquareSpaceGroupBasis) = PXPSquareSpaceGroupBasis(
    copy(b.dgt),
    b.I,
    b.R,
    b.B,
    b.n,
    b.point_group,
    b.constrained_count,
    b.perms,
    b.perm_tables,
    b.lookup,
    b.mask,
    b.row_mask,
)

EDKit.order(b::PXPSquareSpaceGroupBasis) = length(b.perms)

"""
    PXPEEDBenchmarkConfig(n; kwargs...)

Configuration for a finite periodic square-lattice PXP ED benchmark.

Keyword defaults are chosen for a short smoke trajectory:
`total_time = 0.1`, `dt = 0.01`, `measure_every = 1`, `initial_state = :down`,
`point_group = true`, `use_sparse = true`, `tol = 1e-10`, `m_init = 30`,
`m_max = 60`, and `extend_step = 10`.
"""
struct PXPEEDBenchmarkConfig
    n::Int
    total_time::Float64
    dt::Float64
    measure_every::Int
    initial_state::Symbol
    point_group::Bool
    use_sparse::Bool
    tol::Float64
    m_init::Int
    m_max::Int
    extend_step::Int
end

function PXPEEDBenchmarkConfig(
    n::Integer;
    total_time::Real = 0.1,
    dt::Real = 0.01,
    measure_every::Integer = 1,
    initial_state::Symbol = :down,
    point_group::Bool = true,
    use_sparse::Bool = true,
    tol::Real = 1e-10,
    m_init::Integer = 30,
    m_max::Integer = 60,
    extend_step::Integer = 10,
)
    n_int = Int(n)
    total = Float64(total_time)
    step = Float64(dt)
    cadence = Int(measure_every)
    n_int >= 3 || throw(ArgumentError("n must be at least 3 for periodic square PXP"))
    n_int^2 <= _MAX_UINT64_SITES ||
        throw(ArgumentError("PXP ED basis currently supports at most $_MAX_UINT64_SITES sites"))
    isfinite(total) && total >= 0 ||
        throw(ArgumentError("total_time must be finite and nonnegative"))
    isfinite(step) && step > 0 || throw(ArgumentError("dt must be finite and positive"))
    _is_integer_multiple(total, step) ||
        throw(ArgumentError("total_time must be an integer multiple of dt"))
    cadence >= 1 || throw(ArgumentError("measure_every must be at least 1"))
    tol_f = Float64(tol)
    isfinite(tol_f) && tol_f > 0 || throw(ArgumentError("tol must be finite and positive"))
    m_init_i = Int(m_init)
    m_max_i = Int(m_max)
    extend_i = Int(extend_step)
    m_init_i >= 1 || throw(ArgumentError("m_init must be at least 1"))
    m_max_i >= m_init_i || throw(ArgumentError("m_max must be at least m_init"))
    extend_i >= 1 || throw(ArgumentError("extend_step must be at least 1"))
    initial_state in (:down, :all_down) ||
        throw(ArgumentError("supported PXP ED initial states are :down and :all_down"))
    return PXPEEDBenchmarkConfig(
        n_int,
        total,
        step,
        cadence,
        initial_state,
        point_group,
        use_sparse,
        tol_f,
        m_init_i,
        m_max_i,
        extend_i,
    )
end

"""
    PXPEEDSample

One measured time sample from a finite PXP ED benchmark trajectory.
"""
struct PXPEEDSample
    step::Int
    time::Float64
    norm::Float64
    return_probability::Float64
    excitation_density::Float64
end

"""
    PXPEEDBenchmarkResult

Result record for `run_pxp_ed_benchmark`, including basis size, Hamiltonian
storage information, trajectory samples, and EDKit Krylov diagnostics.
"""
struct PXPEEDBenchmarkResult
    lattice_size::Tuple{Int,Int}
    basis_dimension::Int
    constrained_dimension::Int
    group_order::Int
    point_group::Bool
    hamiltonian_nnz::Union{Nothing,Int}
    samples::Vector{PXPEEDSample}
    diagnostics::EDKit.KrylovEvolutionDiagnostics
end

function _is_integer_multiple(total::Float64, dt::Float64)
    total == 0.0 && return true
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

function _site_index(n::Integer, x::Integer, y::Integer)
    return x + n * y + 1
end

function _row_masks(n::Integer)
    last_bit = 1 << (n - 1)
    limit = 1 << n
    masks = Int[]
    for mask = 0:(limit - 1)
        has_linear_adjacent = (mask & (mask << 1)) != 0
        has_wrap_adjacent = (mask & 1) != 0 && (mask & last_bit) != 0
        (has_linear_adjacent || has_wrap_adjacent) || push!(masks, mask)
    end
    return masks
end

function _row_occ_table(n::Integer)
    nsites = n^2
    table = zeros(UInt64, n, 1 << n)
    for y = 0:(n - 1), row = 0:((1 << n) - 1)
        occ = UInt64(0)
        for x = 0:(n - 1)
            if !iszero(row & (1 << x))
                site = _site_index(n, x, y)
                occ |= UInt64(1) << (nsites - site)
            end
        end
        table[y + 1, row + 1] = occ
    end
    return table
end

function _space_group_perms(n::Integer, point_group::Bool)
    point_ops = if point_group
        (
            (x, y) -> (x, y),
            (x, y) -> (y, n - 1 - x),
            (x, y) -> (n - 1 - x, n - 1 - y),
            (x, y) -> (n - 1 - y, x),
            (x, y) -> (n - 1 - x, y),
            (x, y) -> (x, n - 1 - y),
            (x, y) -> (y, x),
            (x, y) -> (n - 1 - y, n - 1 - x),
        )
    else
        ((x, y) -> (x, y),)
    end

    perms = Vector{Int}[]
    seen = Set{Tuple{Vararg{Int}}}()
    for op in point_ops, dy = 0:(n - 1), dx = 0:(n - 1)
        perm = Vector{Int}(undef, n^2)
        for y = 0:(n - 1), x = 0:(n - 1)
            x2, y2 = op(x, y)
            perm[_site_index(n, x, y)] = _site_index(n, mod(x2 + dx, n), mod(y2 + dy, n))
        end
        key = Tuple(perm)
        if !(key in seen)
            push!(seen, key)
            push!(perms, perm)
        end
    end
    return perms
end

function _permutation_tables(n::Integer, perms::Vector{Vector{Int}})
    nsites = n^2
    tables = Matrix{UInt64}[]
    for perm in perms
        table = zeros(UInt64, n, 1 << n)
        for y = 0:(n - 1), chunk = 0:((1 << n) - 1)
            out = UInt64(0)
            for x = 0:(n - 1)
                if !iszero(chunk & (1 << (n - 1 - x)))
                    src = _site_index(n, x, y)
                    dst = perm[src]
                    out |= UInt64(1) << (nsites - dst)
                end
            end
            table[y + 1, chunk + 1] = out
        end
        push!(tables, table)
    end
    return tables
end

@inline function _apply_permutation(bits::UInt64, table::Matrix{UInt64}, n::Int, row_mask::UInt64)
    out = UInt64(0)
    @inbounds for y = 0:(n - 1)
        shift = n * (n - y - 1)
        chunk = Int((bits >> shift) & row_mask)
        out |= table[y + 1, chunk + 1]
    end
    return out
end

function _canonical_digit_bits(bits::UInt64, tables::Vector{Matrix{UInt64}}, n::Int, row_mask::UInt64)
    best = bits
    @inbounds for table in tables
        transformed = _apply_permutation(bits, table, n, row_mask)
        transformed < best && (best = transformed)
    end
    return best
end

function _canonical_representative_data(
    bits::UInt64,
    tables::Vector{Matrix{UInt64}},
    n::Int,
    row_mask::UInt64,
)
    stabilizer = 0
    @inbounds for table in tables
        transformed = _apply_permutation(bits, table, n, row_mask)
        transformed < bits && return false, 0
        transformed == bits && (stabilizer += 1)
    end
    return true, stabilizer
end

function _digits_to_bits(dgt::AbstractVector{<:Integer})
    bits = UInt64(0)
    for digit in dgt
        digit == 0 || digit == 1 || throw(ArgumentError("PXP ED basis digits must be 0 or 1"))
        bits = (bits << 1) | UInt64(digit)
    end
    return bits
end

function _is_constrained_digits(dgt::AbstractVector{<:Integer}, n::Integer)
    length(dgt) == n^2 || throw(ArgumentError("digit vector length must be n^2"))
    for y = 0:(n - 1), x = 0:(n - 1)
        site = _site_index(n, x, y)
        if dgt[site] == 0
            right = _site_index(n, mod(x + 1, n), y)
            up = _site_index(n, x, mod(y + 1, n))
            (dgt[right] == 0 || dgt[up] == 0) && return false
        end
    end
    return true
end

function _enumerate_pbc_independent_sets(n::Integer, visit)
    masks = _row_masks(n)
    row_occ = _row_occ_table(n)
    compatible = Dict(mask => Int[] for mask in masks)
    for prev in masks, next in masks
        iszero(prev & next) && push!(compatible[prev], next)
    end

    rows = zeros(Int, n)
    function rec(y::Int)
        if y > n
            iszero(rows[end] & rows[1]) || return nothing
            occ = UInt64(0)
            @inbounds for row_y = 1:n
                occ |= row_occ[row_y, rows[row_y] + 1]
            end
            visit(occ)
            return nothing
        end

        candidates = y == 1 ? masks : compatible[rows[y - 1]]
        for mask in candidates
            rows[y] = mask
            rec(y + 1)
        end
        return nothing
    end
    rec(1)
    return nothing
end

"""
    pxp_ed_space_group_basis(n; point_group = true)

Construct the fully symmetric EDKit-compatible PXP basis for an `n x n`
periodic square lattice.

With `point_group = true`, the basis keeps one representative per orbit under
translations and the square point group. With `point_group = false`, only
translations are used. The implementation enumerates periodic hard-square
states directly instead of scanning the full `2^(n^2)` product space.
"""
function pxp_ed_space_group_basis(n::Integer; point_group::Bool = true)
    n_int = Int(n)
    n_int >= 3 || throw(ArgumentError("n must be at least 3 for periodic square PXP"))
    nsites = n_int^2
    nsites <= _MAX_UINT64_SITES ||
        throw(ArgumentError("PXP ED basis currently supports at most $_MAX_UINT64_SITES sites"))

    perms = _space_group_perms(n_int, point_group)
    tables = _permutation_tables(n_int, perms)
    mask = (UInt64(1) << nsites) - UInt64(1)
    row_mask = (UInt64(1) << n_int) - UInt64(1)
    entries = Tuple{Int,Float64}[]
    constrained_count = 0

    _enumerate_pbc_independent_sets(n_int, occ_bits -> begin
        constrained_count += 1
        digit_bits = mask ⊻ occ_bits
        keep, stabilizer = _canonical_representative_data(digit_bits, tables, n_int, row_mask)
        if keep
            push!(entries, (Int(digit_bits + UInt64(1)), sqrt(length(perms) * stabilizer)))
        end
        nothing
    end)

    sort!(entries; by = first)
    reps = [entry[1] for entry in entries]
    norms = [entry[2] for entry in entries]
    lookup = Dict{Int,Int}(rep => idx for (idx, rep) in pairs(reps))
    return PXPSquareSpaceGroupBasis(
        zeros(Int, nsites),
        reps,
        norms,
        2,
        n_int,
        point_group,
        constrained_count,
        perms,
        tables,
        lookup,
        mask,
        row_mask,
    )
end

"""
    pxp_ed_constrained_count(basis)

Return the number of PBC blockade-allowed product states before symmetry
reduction for a PXP ED basis.
"""
pxp_ed_constrained_count(basis::PXPSquareSpaceGroupBasis) = basis.constrained_count

"""
    pxp_ed_group_order(basis)

Return the finite symmetry-group order used by a PXP ED basis.
"""
pxp_ed_group_order(basis::PXPSquareSpaceGroupBasis) = length(basis.perms)

"""
    pxp_ed_boundary_condition(basis)

Return the boundary condition represented by a PXP ED basis.
"""
pxp_ed_boundary_condition(::PXPSquareSpaceGroupBasis) = :periodic

"""
    pxp_ed_symmetry_sector(basis)

Return the symmetry sector represented by a PXP ED basis.
"""
pxp_ed_symmetry_sector(basis::PXPSquareSpaceGroupBasis) =
    basis.point_group ? :fully_symmetric_space_group : :translation_symmetric

"""
    pxp_ed_observable_scope(basis)

Return the observable scope supported by the current PBC symmetry-reduced ED
basis. The value is global because local and central-region observables do not
preserve the selected symmetric sector.
"""
pxp_ed_observable_scope(::PXPSquareSpaceGroupBasis) = :pbc_global_site_average

"""
    pxp_ed_reference_label(basis)

Return a stable machine-readable label for the ED reference observable.
"""
pxp_ed_reference_label(::PXPSquareSpaceGroupBasis) = "finite_pbc_global_density"

"""
    pxp_ed_site_density_operator(basis, site)

Construct a site-density operator when the supplied basis supports local
observables. The current symmetry-reduced PBC basis rejects this request because
it would be a projected group average, not a site observable.
"""
function pxp_ed_site_density_operator(::PXPSquareSpaceGroupBasis, site::Integer)
    site >= 1 || throw(ArgumentError("site must be positive"))
    throw(ArgumentError("site density is not available in the symmetry-reduced PBC ED basis"))
end

"""
    pxp_ed_region_density_operator(basis, sites)

Construct a region-density operator when the supplied basis supports local
regions. The current symmetry-reduced PBC basis rejects this request because
there is no central region in a fully symmetric periodic basis.
"""
function pxp_ed_region_density_operator(::PXPSquareSpaceGroupBasis, sites)
    isempty(collect(sites)) && throw(ArgumentError("region sites must be nonempty"))
    throw(ArgumentError("region density is not available in the symmetry-reduced PBC ED basis"))
end

function EDKit.index(basis::PXPSquareSpaceGroupBasis, dgt::AbstractVector)
    _is_constrained_digits(dgt, basis.n) || return 0.0, 1
    bits = _digits_to_bits(dgt)
    canonical = _canonical_digit_bits(bits, basis.perm_tables, basis.n, basis.row_mask)
    pos = get(basis.lookup, Int(canonical + UInt64(1)), 0)
    pos == 0 && return 0.0, 1
    return basis.R[pos], pos
end

"""
    pxp_ed_initial_state(basis; state = :down)

Return a normalized coordinate vector for a symmetry-compatible benchmark
initial state. Currently supported states are `:down` and `:all_down`.
"""
function pxp_ed_initial_state(basis::PXPSquareSpaceGroupBasis; state::Symbol = :down)
    state === :down || state === :all_down ||
        throw(ArgumentError("supported PXP ED initial states are :down and :all_down"))
    dgt = ones(Int, length(basis.dgt))
    coeff, pos = EDKit.index(basis, dgt)
    iszero(coeff) && throw(ArgumentError("all-down state is not present in the supplied basis"))
    psi = zeros(ComplexF64, size(basis, 1))
    psi[pos] = 1
    return psi
end

"""
    pxp_ed_hamiltonian_operator(basis)

Construct the EDKit `Operator` for the constrained square-lattice PXP
Hamiltonian `sum_i X_i` in the supplied PBC symmetry-reduced basis.
"""
function pxp_ed_hamiltonian_operator(basis::PXPSquareSpaceGroupBasis)
    nsites = length(basis.dgt)
    return EDKit.operator(fill(EDKit.spin("X"), nsites), collect(1:nsites), basis)
end

"""
    sparse_pxp_ed_hamiltonian(basis)

Materialize the PXP ED Hamiltonian as a sparse matrix. This is usually faster
for repeated Krylov matvecs, but it costs more memory than the matrix-free
EDKit `Operator`.
"""
sparse_pxp_ed_hamiltonian(basis::PXPSquareSpaceGroupBasis) =
    sparse(pxp_ed_hamiltonian_operator(basis))

function _pxp_ed_density_operator(basis::PXPSquareSpaceGroupBasis)
    nsites = length(basis.dgt)
    n_up = spdiagm(0 => [1.0, 0.0])
    return (1 / nsites) * EDKit.operator(fill(n_up, nsites), collect(1:nsites), basis)
end

function _sample(step::Integer, time::Real, psi, psi0, density_op)
    norm_value = norm(psi)
    return PXPEEDSample(
        Int(step),
        Float64(time),
        norm_value,
        abs2(dot(psi0, psi)),
        real(dot(psi, density_op * psi)),
    )
end

"""
    run_pxp_ed_benchmark(config)

Run a short-time finite PBC square-lattice PXP ED benchmark with EDKit adaptive
Krylov evolution. The Hamiltonian is built from an EDKit `Operator`; by default
it is materialized as a sparse matrix before time evolution.
"""
function run_pxp_ed_benchmark(config::PXPEEDBenchmarkConfig)
    basis = pxp_ed_space_group_basis(config.n; point_group = config.point_group)
    h_operator = pxp_ed_hamiltonian_operator(basis)
    hamiltonian = config.use_sparse ? sparse(h_operator) : h_operator
    hamiltonian_nnz = hamiltonian isa SparseMatrixCSC ? nnz(hamiltonian) : nothing
    density_op = sparse(_pxp_ed_density_operator(basis))
    psi0 = pxp_ed_initial_state(basis; state = config.initial_state)

    cache = EDKit.KrylovEvolutionCache(
        hamiltonian,
        psi0;
        tol = config.tol,
        m_init = config.m_init,
        m_max = config.m_max,
        extend_step = config.extend_step,
    )

    samples = PXPEEDSample[_sample(0, 0.0, psi0, psi0, density_op)]
    nsteps = round(Int, config.total_time / config.dt)
    for step = 1:nsteps
        if step % config.measure_every == 0 || step == nsteps
            time = step * config.dt
            psi = EDKit.timeevolve!(cache, time)
            push!(samples, _sample(step, time, psi, psi0, density_op))
        end
    end

    return PXPEEDBenchmarkResult(
        (config.n, config.n),
        size(basis, 1),
        pxp_ed_constrained_count(basis),
        pxp_ed_group_order(basis),
        config.point_group,
        hamiltonian_nnz,
        samples,
        cache.diagnostics,
    )
end

function _sample_data(sample::PXPEEDSample)
    return (;
        step = sample.step,
        time = sample.time,
        norm = sample.norm,
        return_probability = sample.return_probability,
        excitation_density = sample.excitation_density,
    )
end

function _diagnostics_data(diagnostics::EDKit.KrylovEvolutionDiagnostics)
    return (;
        basis_builds = diagnostics.basis_builds,
        basis_extensions = diagnostics.basis_extensions,
        restarts = diagnostics.restarts,
        matvecs = diagnostics.matvecs,
        total_times_served = diagnostics.total_times_served,
        max_dim_used = diagnostics.max_dim_used,
        accepted_intervals = diagnostics.accepted_intervals,
    )
end

"""
    write_pxp_ed_benchmark_json(result, path)

Write a `PXPEEDBenchmarkResult` as a JSON record.
"""
function write_pxp_ed_benchmark_json(result::PXPEEDBenchmarkResult, path::AbstractString)
    data = (;
        lattice_size = result.lattice_size,
        basis_dimension = result.basis_dimension,
        constrained_dimension = result.constrained_dimension,
        group_order = result.group_order,
        point_group = result.point_group,
        hamiltonian_nnz = result.hamiltonian_nnz,
        samples = _sample_data.(result.samples),
        diagnostics = _diagnostics_data(result.diagnostics),
    )
    open(path, "w") do io
        JSON3.write(io, data)
        write(io, '\n')
    end
    return path
end

end
