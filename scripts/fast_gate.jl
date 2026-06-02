#!/usr/bin/env julia
# Fast physics-correctness gate — the per-iteration heartbeat for the
# review -> test -> benchmark -> improve loop.
#
# Runs only the tests that would catch a Stage-1/Stage-2 physics regression,
# in ONE process (no Distributed), so the ~75s first-call compile tax is paid
# once rather than per worker. Everything else (Aqua, schema/serialization
# smoke, real-CTM, TFIM references) belongs in the nightly full suite.
#
# Usage:
#   julia --project=test scripts/fast_gate.jl
#   julia --project=test scripts/fast_gate.jl test_star_simple_update.jl   # subset
#
# Exit code is nonzero if any selected file fails, so it is safe as a CI/loop gate.

using Test
using Printf
using SquarePXPDynamics

const FAST_FILES = [
    "test_star_simple_update.jl",   # D=1 dense square-star reference exactness
    "test_ipeps_evolution.jl",      # Trotter schedule + normalization ledger
    "test_pxp_d2_localization.jl",  # D=2 exact-finite vs simple separation
    "test_pxp_validation.jl",       # D=1 iPEPS-vs-ED density match + report shape
    "test_ipeps_compression.jl",    # bond compression / truncation bookkeeping
]

const TEST_DIR = joinpath(dirname(@__DIR__), "test")
const SELECTED = isempty(ARGS) ? FAST_FILES : ARGS

failed = String[]
total = @elapsed for f in SELECTED
    print(stderr, "▶ fast-gate: $f ... ")
    ok = true
    t = @elapsed try
        Base.include(Main, joinpath(TEST_DIR, f))
    catch e
        ok = false
        push!(failed, f)
        showerror(stderr, e, catch_backtrace())
        println(stderr)
    end
    @printf(stderr, "%s  %.1fs\n", ok ? "PASS" : "FAIL", t)
end

println(stderr, "="^56)
@printf(stderr, "fast-gate %s  (%d files, %.1fs total)\n",
        isempty(failed) ? "PASS ✓" : "FAIL ✗", length(SELECTED), total)
if !isempty(failed)
    println(stderr, "  failed: ", join(failed, ", "))
    exit(1)
end
