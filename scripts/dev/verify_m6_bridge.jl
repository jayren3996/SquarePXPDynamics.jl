using SquarePXPDynamics
using SquarePXPDynamics: PeriodicSquareUnitCell, SquareCoord
using SquarePXPDynamics.Observables: measure_simple
using SquarePXPDynamics.IPEPSEvolution: TrotterParams, evolve!
using Printf

approx(a, b; atol=1e-9) = isapprox(a, b; atol=atol)

function cmp_summaries(name, a, b; atol=1e-9)
    fields = (:density, :density_even, :density_odd, :blockade_violation,
              :pxp_energy_density, :mean_bond_entropy, :max_bond_entropy)
    ok = true
    for f in fields
        va = getfield(a, f); vb = getfield(b, f)
        d = abs(va - vb)
        flag = d <= atol ? "ok" : "MISMATCH"
        d <= atol || (ok = false)
        @printf("  %-22s native=% .10f legacy=% .10f  |Δ|=%.2e  %s\n", f, va, vb, d, flag)
    end
    println("=> $name: ", ok ? "PASS" : "FAIL")
    return ok
end

println("== D=1 product (:down) parity ==")
cell = PeriodicSquareUnitCell(4, 4)
nat = product_pxp_pepskit_state(cell; state=:down, D=1)
leg = product_square_ipeps(cell; state=:down, maxdim=1)
ok1 = cmp_summaries("D=1 down", measure_simple(nat), measure_simple(leg))

println("== D=1 checkerboard parity ==")
natc = checkerboard_pxp_pepskit_state(cell; excited_on=:even, D=1)
legc = checkerboard_square_ipeps(cell; excited_on=:even, maxdim=1)
ok2 = cmp_summaries("D=1 checkerboard", measure_simple(natc), measure_simple(legc))

println("== D>1 converter parity (evolve checkerboard legacy, convert, compare) ==")
# Evolve a checkerboard legacy state to grow D>1 with NONZERO density/blockade/
# energy, then convert to native and compare simple observables in the SAME
# Vidal gauge. This exercises the squared-Schmidt cut-leg weighting on the
# single-, two-, and five-site patches.
cell5 = PeriodicSquareUnitCell(5, 5)
psi = product_square_ipeps(cell5; state=:down, maxdim=4)
tp = TrotterParams(0.3, 2, :real, 4, 1e-12; schedule=:five_color)
evolve!(psi, 0.6; params=tp, projected=false)   # non-projected → nonzero blockade
maxD = maximum(length(SquarePXPDynamics.link_weight(psi, c, :right)) for c in cell5.reps)
sl = measure_simple(psi)
println("  grown legacy max right-bond D = ", maxD,
        "  (legacy density=", round(sl.density; digits=4),
        " blockade=", round(sl.blockade_violation; digits=4),
        " energy=", round(sl.pxp_energy_density; digits=4), ")")
conv = pxpipeps_from_square_ipeps(psi)
ok3 = cmp_summaries("D>1 converted", measure_simple(conv), sl; atol=1e-8)

println("== native evolve! -> EvolutionLog ==")
st = product_pxp_pepskit_state(cell5; state=:down, D=1)
log = evolve!(st, 0.1; params=TrotterParams(0.05, 2, :real, 4, 1e-12; schedule=:five_color))
@printf("  nsteps=%d max_truncerr=%.3e max_bond_entropy=%.3e log_norm_delta=%.3e\n",
        log.nsteps, log.max_truncerr, log.max_bond_entropy, log.log_norm_delta)
ok4 = (log.nsteps == 2) && isfinite(log.max_truncerr) && isfinite(log.max_bond_entropy)
println("=> native evolve!: ", ok4 ? "PASS" : "FAIL")

println("\nALL: ", all((ok1, ok2, ok3, ok4)) ? "PASS" : "FAIL")
