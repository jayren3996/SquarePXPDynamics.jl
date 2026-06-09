using SquarePXPDynamics
using Printf

# Lean diagnostic: D=2 only, 4x4 Neel. At each 0.2-sample record exact-oracle density
# AND the per-sample worst simple-update truncation error (from EvolutionLog.max_truncerr).
# Question: during the revival RISE (t~1.4-2.6) where the lag lives, is the truncation
# discarding meaningful weight (=> a metric could matter) or ~zero (=> missing correlation)?
TS = collect(0.0:0.2:2.8)
ED = [0.5, 0.48026525133620573, 0.4241792325244558, 0.34070816732757286,
      0.2442464623167832, 0.15549264928630516, 0.10057443671942552,
      0.10056019790056131, 0.1555774954125962, 0.24370069054283897,
      0.3373453880372616, 0.4159980265031265, 0.46695735678274974,
      0.4830586210094464, 0.46229670220636837]
DT=0.02; SPS=10
D=2
psi = checkerboard_square_ipeps(PeriodicSquareUnitCell(4,4); excited_on=:even, maxdim=D)
params = TrotterParams(DT, 2, :real, D, 1e-12; schedule=:serial)
@printf("%-6s %-10s %-10s %-12s\n","t","n_iPEPS","|err|","sample_truncerr")
n0 = exact_density_finite(psi; max_sites=16)
@printf("%-6.1f %-10.5f %-10.4e %-12.4e\n", 0.0, n0, abs(n0-ED[1]), 0.0)
for i in 2:length(TS)
    log = evolve!(psi, 0.2; params=params)   # one 0.2 chunk; max_truncerr over its 10 steps
    n = exact_density_finite(psi; max_sites=16)
    @printf("%-6.1f %-10.5f %-10.4e %-12.4e\n", TS[i], n, abs(n-ED[i]), log.max_truncerr)
    flush(stdout)
end
println("DONE")
