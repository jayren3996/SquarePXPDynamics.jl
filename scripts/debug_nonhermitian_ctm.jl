#!/usr/bin/env julia
#
# Minimal reproduction of the non-Hermitian CTM observable seen in the v1
# audit at iterations>1. Compares two states that should be physically
# equivalent:
#   A. Single evolve!(psi, T)
#   B. N × evolve!(psi, T/N) (the loop pattern scarfinder! uses)
# For each, prints the raw complex CTM density and PXP energy expectation
# values, with a softer atol so non-Hermitian artifacts surface as numbers
# rather than exceptions.

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using ITensors: scalar
using PEPSKit: expectation_value
using Printf
using SquarePXPDynamics
import SquarePXPDynamics: PEPSKitMeasurements as PM

configure_ctm_threading_from_env!()

function _raw_density_expectation(psi)
    ctx = pepskit_ctmrg_context(psi; params = default_ctmrg_params(; chi = 16))
    op = PM._cached_density_sum_operator(psi, ctx; sublattice = nothing)
    raw = expectation_value(ctx.peps, op, ctx.env)
    reps = PM._sublattice_reps(psi.unitcell, nothing)
    return raw / length(reps)
end

function _raw_pxp_energy_expectation(psi)
    ctx = pepskit_ctmrg_context(psi; params = default_ctmrg_params(; chi = 16))
    op = PM._cached_pxp_energy_operator(psi, ctx)
    raw = expectation_value(ctx.peps, op, ctx.env)
    reps = PM._sublattice_reps(psi.unitcell, nothing)
    return raw / length(reps)
end

function _build_single(total_time, dt, D)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    trotter = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
    evolve!(psi, total_time; params = trotter)
    return psi
end

function _build_iterated(total_time, dt, D, N)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    trotter = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
    step = total_time / N
    for _ in 1:N
        evolve!(psi, step; params = trotter)
    end
    return psi
end

function _build_iterated_with_intermediate_ctm(total_time, dt, D, N)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    trotter = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
    step = total_time / N
    for _ in 1:N
        evolve!(psi, step; params = trotter)
        # mimic scarfinder!: measure between iterations
        measure_simple(psi)
        measure_ctm(psi; params = default_ctmrg_params(; chi = 16))
    end
    return psi
end

function _build_iterated_with_simple(total_time, dt, D, N)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    trotter = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
    step = total_time / N
    for _ in 1:N
        evolve!(psi, step; params = trotter)
        measure_simple(psi)
    end
    return psi
end

function _report(label, value)
    z = ComplexF64(value)
    @printf("  %-40s  Re = %+.9e   Im = %+.3e\n", label, real(z), imag(z))
end

total_time = 0.48  # 24 dt-steps at dt=0.02; divisible by N in {2,3,4,6,8,12}
dt = 0.02
D = 2

println("Single evolve!(psi, $total_time):")
psi_a = _build_single(total_time, dt, D)
_report("density", _raw_density_expectation(psi_a))
_report("pxp_energy_density", _raw_pxp_energy_expectation(psi_a))
println("  log_norm = ", log_norm(psi_a))

for N in (2, 3, 6)
    println("\n$N × evolve!(psi, $(total_time/N)):")
    psi_b = _build_iterated(total_time, dt, D, N)
    _report("density", _raw_density_expectation(psi_b))
    _report("pxp_energy_density", _raw_pxp_energy_expectation(psi_b))
    println("  log_norm = ", log_norm(psi_b))
end

for N in (3, 6)
    println("\n$N × (evolve! + measure_simple) (psi, $(total_time/N) each):")
    psi_c = _build_iterated_with_simple(total_time, dt, D, N)
    _report("density", _raw_density_expectation(psi_c))
    _report("pxp_energy_density", _raw_pxp_energy_expectation(psi_c))
    println("  log_norm = ", log_norm(psi_c))
end

for N in (3, 6)
    println("\n$N × (evolve! + measure_simple + measure_ctm) (psi, $(total_time/N) each):")
    psi_d = _build_iterated_with_intermediate_ctm(total_time, dt, D, N)
    _report("density", _raw_density_expectation(psi_d))
    _report("pxp_energy_density", _raw_pxp_energy_expectation(psi_d))
    println("  log_norm = ", log_norm(psi_d))
end

function _build_via_scarfinder(total_time, dt, D, N, chi)
    cell = PeriodicSquareUnitCell(3, 3)
    psi = product_square_ipeps(cell; state = :down, maxdim = D)
    trotter = TrotterParams(dt, 1, :real, D, 1e-12; schedule = :serial)
    sf_params = ScarFinderParams(total_time, trotter, N, Inf, Inf, Inf, false)
    ctm_window = (
        PEPSKitCTMRGParams(max(8, chi - 4), 1e-8, 100, 0; seed = 0),
        PEPSKitCTMRGParams(chi, 1e-8, 100, 0; seed = 0),
    )
    backend = TrustedCTMBackend(ctm_window, CTMTrustPolicy())
    scarfinder!(psi, sf_params; measurement = backend, ctm_every = 1)
    return psi
end

for N in (3, 6)
    println("\nscarfinder!(N=$N, ctm_every=1, chi=[12,16]):")
    psi_e = _build_via_scarfinder(total_time, dt, D, N, 16)
    _report("density", _raw_density_expectation(psi_e))
    _report("pxp_energy_density", _raw_pxp_energy_expectation(psi_e))
    println("  log_norm = ", log_norm(psi_e))
end
