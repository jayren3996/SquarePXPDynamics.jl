module FiniteMPSTFIMReference

using ITensorMPS
using LinearAlgebra

export FiniteMPSTFIMMetadata, FiniteMPSTFIMSample, FiniteMPSTFIMResult
export finite_mps_site_index, finite_mps_square_lattice_bonds
export run_finite_mps_tfim_reference

"""
    FiniteMPSTFIMMetadata

Configuration metadata for a finite open-boundary square-lattice TFIM MPS
reference run.
"""
struct FiniteMPSTFIMMetadata
    Lx::Int
    Ly::Int
    boundary::Symbol
    method::Symbol
    J::Float64
    h::Float64
    initial_state::Symbol
    total_time::Float64
    dt::Float64
    measure_every::Int
    maxdim::Int
    cutoff::Float64
end

"""
    FiniteMPSTFIMSample

One finite-MPS TFIM reference sample with global observables, energy density,
maximum MPS bond dimension, and norm.
"""
struct FiniteMPSTFIMSample
    step::Int
    time::Float64
    mean_x::Float64
    mean_y::Float64
    mean_z::Float64
    zz_right::Float64
    zz_up::Float64
    energy_density::Float64
    maxlinkdim::Int
    norm::Float64
end

"""
    FiniteMPSTFIMResult

Complete finite-MPS TFIM benchmark trajectory returned by
[`run_finite_mps_tfim_reference`](@ref).
"""
struct FiniteMPSTFIMResult
    metadata::FiniteMPSTFIMMetadata
    samples::Vector{FiniteMPSTFIMSample}
end

function _validate_dimensions(Lx::Integer, Ly::Integer)
    nx = Int(Lx)
    ny = Int(Ly)
    nx >= 2 || throw(ArgumentError("Lx must be at least 2"))
    ny >= 2 || throw(ArgumentError("Ly must be at least 2"))
    return nx, ny
end

"""
    finite_mps_site_index(x, y, Lx, Ly)

Return the one-dimensional snake-MPS site index for open square-lattice
coordinate `(x, y)`. Odd rows run left-to-right and even rows right-to-left.
"""
function finite_mps_site_index(x::Integer, y::Integer, Lx::Integer, Ly::Integer)
    nx, ny = _validate_dimensions(Lx, Ly)
    xx = Int(x)
    yy = Int(y)
    1 <= xx <= nx || throw(ArgumentError("x must be in 1:Lx"))
    1 <= yy <= ny || throw(ArgumentError("y must be in 1:Ly"))
    row_offset = (yy - 1) * nx
    row_index = isodd(yy) ? xx : nx - xx + 1
    return row_offset + row_index
end

"""
    finite_mps_square_lattice_bonds(Lx, Ly)

Return open-boundary nearest-neighbor square-lattice bonds as
`(site1, site2, direction)` tuples in the snake-MPS indexing convention.
Directions are `:right` and `:up`.
"""
function finite_mps_square_lattice_bonds(Lx::Integer, Ly::Integer)
    nx, ny = _validate_dimensions(Lx, Ly)
    bonds = Tuple{Int,Int,Symbol}[]
    for y = 1:ny
        for x = 1:nx
            site = finite_mps_site_index(x, y, nx, ny)
            if x < nx
                push!(bonds, (site, finite_mps_site_index(x + 1, y, nx, ny), :right))
            end
            if y < ny
                push!(bonds, (site, finite_mps_site_index(x, y + 1, nx, ny), :up))
            end
        end
    end
    return bonds
end

function _state_name(state::Symbol)
    if state === :z_up || state === :up
        return "0"
    elseif state === :z_down || state === :down
        return "1"
    elseif state === :x_plus
        return "+"
    else
        throw(ArgumentError("initial_state must be :z_up, :z_down, :up, :down, or :x_plus"))
    end
end

function _tfim_opsum(Lx::Int, Ly::Int, J::Float64, h::Float64)
    N = Lx * Ly
    os = OpSum()
    for site = 1:N
        os += -h, "X", site
    end
    for (site1, site2, _) in finite_mps_square_lattice_bonds(Lx, Ly)
        os += -J, "Z", site1, "Z", site2
    end
    return os
end

function _is_integer_multiple(total::Float64, dt::Float64)
    total == 0.0 && return true
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

function _mean_expectation(psi, op::String)
    vals = expect(psi, op)
    return Float64(real(sum(vals) / length(vals)))
end

function _bond_expectation(psi, sites, bonds, dir::Symbol)
    selected = [(i, j) for (i, j, d) in bonds if d === dir]
    isempty(selected) && return 0.0
    os = OpSum()
    for (i, j) in selected
        os += 1.0, "Z", i, "Z", j
    end
    O = MPO(os, sites)
    return Float64(real(inner(psi', O, psi)) / length(selected))
end

function _sample(step::Int, time::Float64, psi, H, sites, bonds)
    energy = inner(psi', H, psi)
    N = length(sites)
    return FiniteMPSTFIMSample(
        step,
        time,
        _mean_expectation(psi, "X"),
        _mean_expectation(psi, "Y"),
        _mean_expectation(psi, "Z"),
        _bond_expectation(psi, sites, bonds, :right),
        _bond_expectation(psi, sites, bonds, :up),
        Float64(real(energy) / N),
        maxlinkdim(psi),
        Float64(real(sqrt(inner(psi, psi)))),
    )
end

"""
    run_finite_mps_tfim_reference(Lx, Ly; kwargs...)

Run an open-boundary `Lx x Ly` square-lattice TFIM reference using a snake-MPS
mapping and ITensorMPS TDVP. The Hamiltonian convention is
`-h * sum_i X_i - J * sum_<ij> Z_i Z_j`.

Keyword defaults are `J = 1.0`, `h = 1.0`, `initial_state = :x_plus`,
`total_time = 0.5`, `dt = 0.05`, `measure_every = 1`, `maxdim = 64`, and
`cutoff = 1e-9`.
"""
function run_finite_mps_tfim_reference(
    Lx::Integer,
    Ly::Integer;
    J::Real = 1.0,
    h::Real = 1.0,
    initial_state::Symbol = :x_plus,
    total_time::Real = 0.5,
    dt::Real = 0.05,
    measure_every::Integer = 1,
    maxdim::Integer = 64,
    cutoff::Real = 1e-9,
    outputlevel::Integer = 0,
)::FiniteMPSTFIMResult
    nx, ny = _validate_dimensions(Lx, Ly)
    coupling = Float64(J)
    field = Float64(h)
    isfinite(coupling) || throw(ArgumentError("J must be finite"))
    isfinite(field) || throw(ArgumentError("h must be finite"))
    total = Float64(total_time)
    step = Float64(dt)
    isfinite(total) && total >= 0 ||
        throw(ArgumentError("total_time must be finite and nonnegative"))
    isfinite(step) && step > 0 || throw(ArgumentError("dt must be finite and positive"))
    _is_integer_multiple(total, step) ||
        throw(ArgumentError("total_time must be an integer multiple of dt"))
    cadence = Int(measure_every)
    cadence >= 1 || throw(ArgumentError("measure_every must be at least 1"))
    dim = Int(maxdim)
    dim >= 1 || throw(ArgumentError("maxdim must be at least 1"))
    trunc_cutoff = Float64(cutoff)
    isfinite(trunc_cutoff) && trunc_cutoff >= 0 ||
        throw(ArgumentError("cutoff must be finite and nonnegative"))

    N = nx * ny
    sites = siteinds("Qubit", N)
    H = MPO(_tfim_opsum(nx, ny, coupling, field), sites)
    psi = MPS(sites, _ -> _state_name(initial_state))
    bonds = finite_mps_square_lattice_bonds(nx, ny)
    nsteps = round(Int, total / step)

    metadata = FiniteMPSTFIMMetadata(
        nx,
        ny,
        :open,
        :tdvp,
        coupling,
        field,
        initial_state,
        total,
        step,
        cadence,
        dim,
        trunc_cutoff,
    )
    samples = FiniteMPSTFIMSample[_sample(0, 0.0, psi, H, sites, bonds)]

    for idx = 1:nsteps
        psi = tdvp(
            H,
            -im * step,
            psi;
            time_step = -im * step,
            cutoff = trunc_cutoff,
            maxdim = dim,
            nsite = 2,
            outputlevel = Int(outputlevel),
        )
        if idx % cadence == 0 || idx == nsteps
            push!(samples, _sample(idx, idx * step, psi, H, sites, bonds))
        end
    end

    return FiniteMPSTFIMResult(metadata, samples)
end

end
