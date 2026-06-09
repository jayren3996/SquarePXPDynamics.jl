using SquarePXPDynamics
using Printf
# D=4 trajectory + truncerr: confirm the lag-vs-truncerr DECOUPLING at the validated-best D
# and quantify the D-scaling of the rise lag (vs the D=2 probe). D=4 dense oracle ~feasible.
TS = collect(0.0:0.2:2.8)
ED = [0.5, 0.48026525133620573, 0.4241792325244558, 0.34070816732757286,
      0.2442464623167832, 0.15549264928630516, 0.10057443671942552,
      0.10056019790056131, 0.1555774954125962, 0.24370069054283897,
      0.3373453880372616, 0.4159980265031265, 0.46695735678274974,
      0.4830586210094464, 0.46229670220636837]
DT=0.02; D=4
psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4,4); excited_on=:even, maxdim=D)
params = TrotterParams(DT, 2, :real, D, 1e-12; schedule=:serial)
@printf("D=4\n%-6s %-10s %-10s %-12s\n","t","n","|err|","truncerr")
n0=exact_density_finite(psi; max_sites=16); @printf("%-6.1f %-10.5f %-10.4e %-12.4e\n",0.0,n0,abs(n0-ED[1]),0.0)
es=Float64[abs(n0-ED[1])]
for i in 2:length(TS)
    log = evolve!(psi, 0.2; params=params)
    n = exact_density_finite(psi; max_sites=16)
    push!(es, abs(n-ED[i]))
    @printf("%-6.1f %-10.5f %-10.4e %-12.4e\n", TS[i], n, abs(n-ED[i]), log.max_truncerr)
    flush(stdout)
end
@printf("RMS=%.4e MAX=%.4e\n", sqrt(sum(es.^2)/length(es)), maximum(es))
println("DONE")
