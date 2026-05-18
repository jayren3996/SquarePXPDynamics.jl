#!/usr/bin/env julia

using Pkg

project_root = dirname(@__DIR__)
Pkg.activate(project_root; io = devnull)

using JSON3
using SquarePXPDynamics

function _env_value(name::String, default::AbstractString)
    value = get(ENV, name, "")
    return isempty(value) ? String(default) : value
end

function _env_bool(name::String, default::Bool)
    value = lowercase(strip(_env_value(name, string(default))))
    value in ("1", "true", "yes", "on") && return true
    value in ("0", "false", "no", "off") && return false
    throw(ArgumentError("$name must be one of 1,true,yes,on,0,false,no,off"))
end

function _env_int(name::String, default::Int)
    return parse(Int, _env_value(name, string(default)))
end

function _env_float(name::String, default::Float64)
    return parse(Float64, _env_value(name, string(default)))
end

function _env_symbol(name::String, default::Symbol)
    return Symbol(_env_value(name, String(default)))
end

function _is_integer_multiple(total::Float64, dt::Float64)
    nsteps = round(Int, total / dt)
    return isapprox(nsteps * dt, total; atol = 1e-12, rtol = 1e-10)
end

function _csv_value(value::Nothing)
    return ""
end

function _csv_value(value)
    return string(value)
end

function _write_csv(path::AbstractString, rows)
    fields = (
        :step,
        :time,
        :density_simple,
        :blockade_violation_simple,
        :pxp_energy_density_simple,
        :exact_finite_density,
        :max_truncerr,
        :mean_bond_entropy,
        :max_bond_entropy,
        :log_norm,
    )
    open(path, "w") do io
        println(io, join(String.(fields), ","))
        for row in rows
            println(io, join((_csv_value(row[field]) for field in fields), ","))
        end
    end
    return path
end

n = _env_int("SQUAREPXP_IPEPS_N", 3)
total_time = _env_float("SQUAREPXP_IPEPS_TOTAL_TIME", 2.0)
dt = _env_float("SQUAREPXP_IPEPS_DT", 0.02)
measure_every = _env_int("SQUAREPXP_IPEPS_MEASURE_EVERY", 1)
D = _env_int("SQUAREPXP_IPEPS_D", 4)
cutoff = _env_float("SQUAREPXP_IPEPS_CUTOFF", 1e-9)
order = _env_int("SQUAREPXP_IPEPS_ORDER", 1)
schedule = _env_symbol("SQUAREPXP_IPEPS_SCHEDULE", :serial)
initial_state = _env_symbol("SQUAREPXP_IPEPS_INITIAL_STATE", :down)
exact_finite = _env_bool("SQUAREPXP_IPEPS_EXACT_FINITE", true)
exact_finite_max_sites = _env_int("SQUAREPXP_IPEPS_EXACT_FINITE_MAX_SITES", 12)

_is_integer_multiple(total_time, dt) ||
    throw(ArgumentError("SQUAREPXP_IPEPS_TOTAL_TIME must be an integer multiple of SQUAREPXP_IPEPS_DT"))
measure_every >= 1 || throw(ArgumentError("SQUAREPXP_IPEPS_MEASURE_EVERY must be at least 1"))

json_out = _env_value(
    "SQUAREPXP_IPEPS_JSON",
    joinpath(project_root, "artifacts", "pxp_ipeps_trajectory.json"),
)
csv_out = _env_value(
    "SQUAREPXP_IPEPS_CSV",
    joinpath(project_root, "artifacts", "pxp_ipeps_trajectory.csv"),
)
mkpath(dirname(json_out))
mkpath(dirname(csv_out))

cell = PeriodicSquareUnitCell(n, n)
psi = product_square_ipeps(cell; state = initial_state, maxdim = D)
trotter = TrotterParams(dt, order, :real, D, cutoff; schedule)
nsteps = round(Int, total_time / dt)

rows = NamedTuple[]

function _sample!(rows, step::Int, time::Float64, psi, evolution)
    push!(
        rows,
        (;
            step,
            time,
            density_simple = density_simple(psi),
            blockade_violation_simple = blockade_violation_simple(psi),
            pxp_energy_density_simple = pxp_energy_density_simple(psi),
            exact_finite_density = exact_finite ?
                exact_density_finite(psi; max_sites = exact_finite_max_sites) : nothing,
            max_truncerr = evolution === nothing ? 0.0 : evolution.max_truncerr,
            mean_bond_entropy = evolution === nothing ? 0.0 : evolution.mean_bond_entropy,
            max_bond_entropy = evolution === nothing ? 0.0 : evolution.max_bond_entropy,
            log_norm = log_norm(psi),
        ),
    )
    return rows
end

_sample!(rows, 0, 0.0, psi, nothing)
for step = 1:nsteps
    evolution = evolve!(psi, dt; params = trotter)
    if step % measure_every == 0 || step == nsteps
        _sample!(rows, step, step * dt, psi, evolution)
    end
end

payload = (;
    config = (;
        n,
        total_time,
        dt,
        measure_every,
        D,
        cutoff,
        order,
        schedule = String(schedule),
        initial_state = String(initial_state),
        exact_finite,
        exact_finite_max_sites,
    ),
    samples = rows,
)

open(json_out, "w") do io
    JSON3.write(io, payload)
    write(io, '\n')
end
_write_csv(csv_out, rows)

println(json_out)
println(csv_out)
