# Probe the actual PEPSKit 0.7 / TensorKit 0.15 API used by the native backend.
# Run from the project env:  julia --project=. scripts/dev/pepskit_api_probe.jl
using PEPSKit
using TensorKit
using InteractiveUtils: subtypes

println("=== PEPSKit version ===")
println(pkgversion(PEPSKit))
println(pkgversion(TensorKit))

names_pepskit = names(PEPSKit)
function has(sym)
    s = Symbol(sym)
    return s in names_pepskit
end

println("\n=== Exported symbols of interest ===")
for s in [:InfinitePEPS, :InfiniteWeightPEPS, :SUWeight, :CTMRGEnv, :leading_boundary,
          :LocalOperator, :expectation_value, :correlator, :correlation_length,
          :PEPSTensor, :SimpleUpdate, :simpleupdate, :absorb_weight, :truncrank,
          :truncerr, :ProductState]
    println(rpad(string(s), 22), " exported=", has(s))
end

println("\n=== SUWeight ===")
println(methods(PEPSKit.SUWeight))
println("fieldnames: ", fieldnames(PEPSKit.SUWeight))

println("\n=== InfiniteWeightPEPS (if present) ===")
if isdefined(PEPSKit, :InfiniteWeightPEPS)
    println(methods(PEPSKit.InfiniteWeightPEPS))
    println("fieldnames: ", fieldnames(PEPSKit.InfiniteWeightPEPS))
end

println("\n=== absorb_weight ===")
if isdefined(PEPSKit, :absorb_weight)
    println(methods(PEPSKit.absorb_weight))
end

println("\n=== InfinitePEPS constructors ===")
println(methods(PEPSKit.InfinitePEPS))

println("\n=== PEPSTensor ===")
if isdefined(PEPSKit, :PEPSTensor)
    println(methods(PEPSKit.PEPSTensor))
end

# Build a tiny D=1 InfinitePEPS and SUWeight to inspect runtime structure
println("\n=== Build tiny 2x2 D=1 state ===")
P = ComplexSpace(2)
V = ComplexSpace(1)
function down_tensor()
    # P <- N E S W  (S,W dualized per existing adapter convention)
    t = zeros(ComplexF64, 2, 1, 1, 1, 1)
    t[2,1,1,1,1] = 1.0  # basis 2 = down
    return TensorMap(t, P ← V ⊗ V ⊗ V' ⊗ V')
end
Lx, Ly = 2, 2
peps = InfinitePEPS([down_tensor() for r in 1:Ly, c in 1:Lx])
println("typeof(peps) = ", typeof(peps))
println("size(peps) ? ")
try; println(size(peps)); catch e; println("size err: ", e); end
println("peps.A size: ", size(peps.A))
println("fieldnames(InfinitePEPS): ", fieldnames(typeof(peps)))

println("\n=== SUWeight construction attempt ===")
# Try the documented constructor signatures
try
    w = SUWeight(peps)
    println("SUWeight(peps) OK: ", typeof(w))
    println("fieldnames: ", fieldnames(typeof(w)))
    try; println("size(w) = ", size(w)); catch e; println(e); end
catch e
    println("SUWeight(peps) failed: ", e)
end

println("\n=== InfiniteWeightPEPS construction attempt ===")
if isdefined(PEPSKit, :InfiniteWeightPEPS)
    try
        wpeps = InfiniteWeightPEPS(peps)
        println("InfiniteWeightPEPS(peps) OK: ", typeof(wpeps))
        println("fieldnames: ", fieldnames(typeof(wpeps)))
    catch e
        println("InfiniteWeightPEPS(peps) failed: ", e)
    end
end

println("\nDONE")
