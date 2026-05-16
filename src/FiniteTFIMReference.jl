module FiniteTFIMReference

using LinearAlgebra

using ..Observables: TFIMObservableSummary
using ..SpinOps: pauli_x, pauli_y, pauli_z, identity2, kron_all, embed_one_site
using ..SquareUnitCells: PeriodicSquareUnitCell, neighbor
using ..StarModels: TFIMStarModel

export FiniteTFIMReferenceSample
export finite_tfim_hamiltonian, finite_tfim_product_state
export measure_finite_tfim, run_finite_tfim_reference

const MAX_DENSE_REFERENCE_SITES = 12

"""
    FiniteTFIMReferenceSample

One exact finite-system TFIM reference sample containing the integer step,
physical time, and convention-matched TFIM observables.
"""
struct FiniteTFIMReferenceSample
    step::Int
    time::Float64
    observables::TFIMObservableSummary
end

function _nsites(cell::PeriodicSquareUnitCell)
    n = cell.Lx * cell.Ly
    n <= MAX_DENSE_REFERENCE_SITES || throw(
        ArgumentError(
            "dense finite TFIM reference supports at most $MAX_DENSE_REFERENCE_SITES sites",
        ),
    )
    return n
end

function _site_index(cell::PeriodicSquareUnitCell, c)
    wrapped = c
    return wrapped.x + (wrapped.y - 1) * cell.Lx
end

function _validate_state_vector(state, nsites::Integer)
    dim = 2^Int(nsites)
    length(state) == dim ||
        throw(ArgumentError("state length must be 2^nsites for the supplied cell"))
    norm_value = norm(state)
    isfinite(norm_value) && norm_value > 0 ||
        throw(ArgumentError("state norm must be finite and nonzero"))
    return Vector{ComplexF64}(state) ./ norm_value
end

function _single_site_state(state::Symbol)
    if state === :z_up || state === :up
        return ComplexF64[1, 0]
    elseif state === :z_down || state === :down
        return ComplexF64[0, 1]
    elseif state === :x_plus
        return ComplexF64[1, 1] ./ sqrt(2)
    else
        throw(ArgumentError("state must be :z_up, :z_down, :up, :down, or :x_plus"))
    end
end

function _expectation(state::Vector{ComplexF64}, op::AbstractMatrix)
    return dot(state, op * state)
end

function _mean_one_site(state::Vector{ComplexF64}, op::AbstractMatrix, nsites::Integer)
    return sum(_expectation(state, embed_one_site(op, site, nsites)) for site = 1:nsites) /
           nsites
end

function _two_site_operator(op1::AbstractMatrix, site1::Integer, op2::AbstractMatrix, site2::Integer, nsites::Integer)
    1 <= site1 <= nsites || throw(ArgumentError("site1 must be in 1:nsites"))
    1 <= site2 <= nsites || throw(ArgumentError("site2 must be in 1:nsites"))
    site1 == site2 && throw(ArgumentError("two-site operator requires distinct sites"))
    return kron_all([
        site == site1 ? op1 : site == site2 ? op2 : identity2() for site = 1:nsites
    ])
end

function _mean_zz(state::Vector{ComplexF64}, cell::PeriodicSquareUnitCell, dir::Symbol)
    dir === :right || dir === :up ||
        throw(ArgumentError("finite TFIM ZZ direction must be :right or :up"))
    nsites = _nsites(cell)
    z = pauli_z()
    total = 0.0 + 0.0im
    for c in cell.reps
        site = _site_index(cell, c)
        nb = _site_index(cell, neighbor(cell, c, dir))
        total += _expectation(state, _two_site_operator(z, site, z, nb, nsites))
    end
    return total / length(cell.reps)
end

function _real_with_imag_abs(value)
    return (real(value), abs(imag(value)))
end

"""
    finite_tfim_hamiltonian(cell, model)

Return the dense finite periodic TFIM Hamiltonian
`-h * sum_i X_i - J * sum_<ij> Z_i Z_j` for right/up square-lattice bonds.
This exact reference is intended for small cells such as `3 x 3`.
"""
function finite_tfim_hamiltonian(
    cell::PeriodicSquareUnitCell,
    model::TFIMStarModel,
)::Matrix{ComplexF64}
    nsites = _nsites(cell)
    dim = 2^nsites
    H = zeros(ComplexF64, dim, dim)
    x = pauli_x()
    z = pauli_z()

    for site = 1:nsites
        H .-= model.h .* embed_one_site(x, site, nsites)
    end

    for c in cell.reps
        site = _site_index(cell, c)
        for dir in (:right, :up)
            nb = _site_index(cell, neighbor(cell, c, dir))
            H .-= model.J .* _two_site_operator(z, site, z, nb, nsites)
        end
    end

    return H
end

"""
    finite_tfim_product_state(cell; state)

Return a normalized finite periodic product-state vector for `cell`.
Supported states are `:z_up`, `:z_down`, and `:x_plus`, with `:up` and
`:down` accepted as aliases.
"""
function finite_tfim_product_state(
    cell::PeriodicSquareUnitCell;
    state::Symbol,
)::Vector{ComplexF64}
    _nsites(cell)
    one_site = _single_site_state(state)
    return kron_all([reshape(one_site, :, 1) for _ in cell.reps])[:, 1]
end

"""
    measure_finite_tfim(state, cell, model)

Measure convention-matched TFIM observables for a normalized finite-state
reference vector. Bond entropy fields are set to zero because the finite exact
reference has no iPEPS truncation diagnostics.
"""
function measure_finite_tfim(
    state,
    cell::PeriodicSquareUnitCell,
    model::TFIMStarModel,
)::TFIMObservableSummary
    nsites = _nsites(cell)
    psi = _validate_state_vector(state, nsites)
    mean_x, x_imag = _real_with_imag_abs(_mean_one_site(psi, pauli_x(), nsites))
    mean_y, y_imag = _real_with_imag_abs(_mean_one_site(psi, pauli_y(), nsites))
    mean_z, z_imag = _real_with_imag_abs(_mean_one_site(psi, pauli_z(), nsites))
    zz_right, zz_right_imag = _real_with_imag_abs(_mean_zz(psi, cell, :right))
    zz_up, zz_up_imag = _real_with_imag_abs(_mean_zz(psi, cell, :up))
    energy = -model.h * mean_x - model.J * (zz_right + zz_up)
    energy_imag = abs(model.h * x_imag) + abs(model.J) * (zz_right_imag + zz_up_imag)
    max_imag = maximum((x_imag, y_imag, z_imag, zz_right_imag, zz_up_imag, energy_imag))

    return TFIMObservableSummary(
        mean_x,
        mean_y,
        mean_z,
        mean_z,
        mean_z,
        zz_right,
        zz_up,
        energy,
        energy,
        0.0,
        x_imag,
        y_imag,
        z_imag,
        max(zz_right_imag, zz_up_imag),
        energy_imag,
        max_imag,
        0.0,
        0.0,
    )
end

function _is_integer_multiple(total::Float64, dt::Float64)
    total == 0.0 && return true
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

"""
    run_finite_tfim_reference(cell, model; initial_state, total_time, dt, measure_every = 1)

Run dense exact finite TFIM time evolution from a product state and return
samples at step 0, every `measure_every` steps, and the final step. This is a
small-cell reference for benchmark validation, not an infinite-system solver.
"""
function run_finite_tfim_reference(
    cell::PeriodicSquareUnitCell,
    model::TFIMStarModel;
    initial_state::Symbol,
    total_time::Real,
    dt::Real,
    measure_every::Integer = 1,
)::Vector{FiniteTFIMReferenceSample}
    total = Float64(total_time)
    step = Float64(dt)
    isfinite(total) && total >= 0 ||
        throw(ArgumentError("total_time must be finite and nonnegative"))
    isfinite(step) && step > 0 || throw(ArgumentError("dt must be finite and positive"))
    _is_integer_multiple(total, step) ||
        throw(ArgumentError("total_time must be an integer multiple of dt"))
    cadence = Int(measure_every)
    cadence >= 1 || throw(ArgumentError("measure_every must be at least 1"))

    H = finite_tfim_hamiltonian(cell, model)
    psi0 = finite_tfim_product_state(cell; state = initial_state)
    nsteps = round(Int, total / step)
    samples = FiniteTFIMReferenceSample[
        FiniteTFIMReferenceSample(0, 0.0, measure_finite_tfim(psi0, cell, model)),
    ]

    for idx = 1:nsteps
        if idx % cadence == 0 || idx == nsteps
            time = idx * step
            psi = exp(-im * time * H) * psi0
            push!(
                samples,
                FiniteTFIMReferenceSample(idx, time, measure_finite_tfim(psi, cell, model)),
            )
        end
    end

    return samples
end

end
