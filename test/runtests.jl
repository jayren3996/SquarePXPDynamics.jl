if abspath(Base.active_project()) == abspath(joinpath(dirname(@__DIR__), "Project.toml"))
    using Pkg
    Pkg.activate(@__DIR__; io = devnull)
end

using Distributed
using Test
using Printf

const TEST_FILES = [
    "test_spinops.jl",
    "test_square_geometry.jl",
    "test_square_pxp.jl",
    "test_star_models.jl",
    "test_square_unitcells.jl",
    "test_square_ipeps.jl",
    "test_square_ipeps_s2.jl",
    "test_star_simple_update.jl",
    "test_ipeps_evolution.jl",
    "test_observables.jl",
    "test_observables_evolved.jl",
    "test_finite_ipeps_observables.jl",
    "test_tfim_schedule_reference.jl",
    "test_pxp_ed_benchmark.jl",
    "test_pxp_larger_d_ed_benchmark.jl",
    "test_ctm_trust.jl",
    "test_pepskit_measurements.jl",
    "test_pxp_validation.jl",
    "test_pxp_d2_localization.jl",
    "test_d_convergence.jl",
    "test_scarfinder.jl",
    "test_candidate_snapshots.jl",
    "test_ipeps_compression.jl",
    "test_scarfinder_audit.jl",
    "test_public_docs.jl",
    "test_aqua.jl",
]

const SELECTED = isempty(ARGS) ? TEST_FILES : ARGS
for f in SELECTED
    f in TEST_FILES || throw(ArgumentError("unknown test file: $f"))
end

const TEST_DIR = @__DIR__
const N_WORKERS = min(
    length(SELECTED),
    parse(Int, get(ENV, "JULIA_TEST_NWORKERS", string(length(SELECTED)))),
)
# Cap BLAS threads per worker so we don't oversubscribe cores
const BLAS_THREADS_PER_WORKER =
    max(1, div(Sys.CPU_THREADS, max(N_WORKERS, 1)))

if N_WORKERS > 1
    addprocs(N_WORKERS; exeflags = `--project=$TEST_DIR`)
end

@everywhere begin
    using Test
    using LinearAlgebra
    using SquarePXPDynamics
    BLAS.set_num_threads($BLAS_THREADS_PER_WORKER)
end

@everywhere const _TEST_DIR = $TEST_DIR

@everywhere function _run_test_file(f)
    ok = true
    err_msg = ""
    elapsed = @elapsed try
        Base.include(Main, joinpath(_TEST_DIR, f))
    catch e
        ok = false
        err_msg = sprint(showerror, e, catch_backtrace())
    end
    return (file = f, ok = ok, err = err_msg, elapsed = elapsed, worker = myid())
end

@info "Running tests in parallel" workers=N_WORKERS files=length(SELECTED) blas_threads_per_worker=BLAS_THREADS_PER_WORKER
results = pmap(_run_test_file, SELECTED)

println("\n", "="^64)
println("Parallel test summary  ($(N_WORKERS) worker(s), $(length(SELECTED)) file(s))")
println("="^64)
for r in sort(results; by = x -> -x.elapsed)
    status = r.ok ? "PASS" : "FAIL"
    @printf("  %-46s %s  %7.2fs  (w%d)\n", r.file, status, r.elapsed, r.worker)
end

failed = filter(r -> !r.ok, results)
if !isempty(failed)
    for r in failed
        println("\n--- FAILED: $(r.file) ---")
        println(r.err)
    end
    error("$(length(failed)) test file(s) failed")
end
