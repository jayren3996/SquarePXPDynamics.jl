#!/usr/bin/env julia

# Exploratory probe for PEPSKit compatibility. This script intentionally uses a
# temporary Julia environment so the package Project.toml/Manifest.toml remain
# untouched.

using Pkg

Pkg.activate(; temp = true)
Pkg.add(["PEPSKit", "TensorKit", "MPSKit", "MatrixAlgebraKit"])

using LinearAlgebra
using PEPSKit
using TensorKit
using MPSKit
using MatrixAlgebraKit

function dependency_version(name)
    for (_, dep) in Pkg.dependencies()
        dep.name == name && return string(dep.version)
    end
    return "unknown"
end

brief(value) = string(typeof(value))
brief(value::Number) = string(value)
brief(value::AbstractArray) = "size=$(size(value)), eltype=$(eltype(value))"
brief(value::NamedTuple) = "NamedTuple$(keys(value))"

function report(f::Function, label)
    print(label, ": ")
    try
        value = f()
        println("ok", value === nothing ? "" : " => $(brief(value))")
        return value
    catch err
        println("failed => ", typeof(err), ": ", sprint(showerror, err))
        return nothing
    end
end

println("Julia ", VERSION)
println("PEPSKit ", dependency_version("PEPSKit"))
println("TensorKit ", dependency_version("TensorKit"))
println("MPSKit ", dependency_version("MPSKit"))

pspace = ComplexSpace(2)
vspace = ComplexSpace(1)

peps = report("instantiate InfinitePEPS 1x1 spin-half D=1") do
    InfinitePEPS(randn, ComplexF64, pspace, vspace; unitcell = (1, 1))
end

report("instantiate InfinitePEPS 2x2 spin-half D=1") do
    InfinitePEPS(randn, ComplexF64, pspace, vspace; unitcell = (2, 2))
end

if peps !== nothing
    report("construct SUWeight") do
        SUWeight(peps)
    end

    report("construct CTMRGEnv with chi=2") do
        CTMRGEnv(randn, ComplexF64, peps, ComplexSpace(2))
    end

    env = report("run minimal CTMRG") do
        env0 = CTMRGEnv(randn, ComplexF64, peps, ComplexSpace(2))
        env, info = leading_boundary(
            env0,
            peps;
            alg = :simultaneous,
            tol = 1.0e-8,
            miniter = 1,
            maxiter = 2,
            verbosity = 0,
            trunc = truncrank(2),
        )
        (env = env, info = info)
    end

    sigma_z = TensorMap(ComplexF64[1 0; 0 -1], pspace ← pspace)
    sigma_x = TensorMap(ComplexF64[0 1; 1 0], pspace ← pspace)
    two_site_xx = sigma_x ⊗ sigma_x
    five_site_dense = reshape(Matrix{ComplexF64}(I, 32, 32), ntuple(_ -> 2, 10))
    five_site_id = TensorMap(five_site_dense, pspace^5 ← pspace^5)
    lattice_1x1 = fill(pspace, 1, 1)
    lattice_3x3 = fill(pspace, 3, 3)
    star_sites = (
        CartesianIndex(2, 2),
        CartesianIndex(2, 3),
        CartesianIndex(1, 2),
        CartesianIndex(2, 1),
        CartesianIndex(3, 2),
    )

    op1 = report("construct 1-site LocalOperator") do
        LocalOperator(lattice_1x1, (CartesianIndex(1, 1),) => sigma_z)
    end

    op2 = report("construct nearest-neighbor 2-site LocalOperator") do
        LocalOperator(lattice_1x1, (CartesianIndex(1, 1), CartesianIndex(1, 2)) => two_site_xx)
    end

    op5 = report("construct custom 5-site star LocalOperator on 3x3 cell") do
        LocalOperator(lattice_3x3, star_sites => five_site_id)
    end

    report("construct custom 5-site LocalCircuit gate") do
        PEPSKit.LocalCircuit(lattice_3x3, star_sites => five_site_id)
    end

    report("trotterize 5-site LocalOperator") do
        PEPSKit.trotterize(op5, 0.01)
    end

    report("time_evolve with 5-site LocalOperator + SimpleUpdate") do
        peps3 = InfinitePEPS(randn, ComplexF64, pspace, vspace; unitcell = (3, 3))
        wts = SUWeight(peps3)
        alg = SimpleUpdate(; trunc = notrunc())
        time_evolve(peps3, op5, 0.01, 1, alg, wts)
    end

    if env !== nothing
        ctm_env = env.env

        report("CTMRG expectation_value 1-site LocalOperator") do
            expectation_value(peps, op1, ctm_env)
        end

        report("CTMRG expectation_value 2-site LocalOperator") do
            expectation_value(peps, op2, ctm_env)
        end

        peps3 = InfinitePEPS(randn, ComplexF64, pspace, vspace; unitcell = (3, 3))
        report("CTMRG expectation_value custom 5-site LocalOperator on 3x3 PEPS") do
            env0 = CTMRGEnv(randn, ComplexF64, peps3, ComplexSpace(2))
            env3, _ = leading_boundary(
                env0,
                peps3;
                alg = :simultaneous,
                tol = 1.0e-8,
                miniter = 1,
                maxiter = 2,
                verbosity = 0,
                trunc = truncrank(2),
            )
            expectation_value(peps3, op5, env3)
        end
    end
end
