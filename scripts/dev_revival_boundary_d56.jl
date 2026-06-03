# D=5,6 revival trajectory via the memory-bounded double-layer boundary oracle.
# The dense path is memory-pathological for entangled D>=5 on a 4x4 torus
# (~2^(N/2)*D^(2Lx); observed ~244 GB), so the boundary contractor is the
# reliable oracle here. Slower per call (one ~D^12 seam SVD per sweep) but bounded.
# Saves after every timepoint; early (collapse) timepoints are cheaper than the
# entangled revival peak.
using SquarePXPDynamics, Printf
const FIO = SquarePXPDynamics.FiniteIPEPSObservables
const TS = collect(0.0:0.2:2.8)
const ED = [0.5,0.48026525133620573,0.4241792325244558,0.34070816732757286,0.2442464623167832,
            0.15549264928630516,0.10057443671942552,0.10056019790056131,0.1555774954125962,
            0.24370069054283897,0.3373453880372616,0.4159980265031265,0.46695735678274974,
            0.4830586210094464,0.46229670220636837]
results = Dict{Int,Vector{Float64}}()
save() = open("/tmp/traj_boundary.json","w") do io
    print(io,"{\"ts\":",TS,",\"ed\":",ED)
    for D in sort(collect(keys(results))); print(io,",\"D",D,"\":",results[D]); end
    print(io,"}")
end
for D in (5,6)
    @printf("=== D=%d ===\n",D); flush(stdout)
    psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4,4); excited_on=:even, maxdim=D)
    params = TrotterParams(0.02,2,:real,D,1e-12; schedule=:serial)
    dens = Float64[]; GC.gc()
    t0=time(); push!(dens, FIO.exact_density_finite(psi; max_sites=16, method=:boundary)); results[D]=copy(dens); save()
    @printf("  t=0.0 n=%.6f (%.0fs)\n", dens[1], time()-t0); flush(stdout)
    for k=2:length(TS)
        for _=1:10; evolve!(psi,0.02; params=params); end
        GC.gc(); t0=time()
        push!(dens, FIO.exact_density_finite(psi; max_sites=16, method=:boundary)); results[D]=copy(dens); save()
        @printf("  t=%.1f n=%.6f |err|=%.4e (%.0fs)\n", TS[k], dens[k], abs(dens[k]-ED[k]), time()-t0); flush(stdout)
    end
    e=abs.(dens.-ED)
    @printf("D=%d rms=%.4e max=%.4e err@2.6=%.4e\n", D, sqrt(sum(e.^2)/15), maximum(e), e[14]); flush(stdout)
end
println("DONE")
