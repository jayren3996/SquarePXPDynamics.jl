using SquarePXPDynamics
const SPD = SquarePXPDynamics

function probe(; n, dt, nsteps, maxdim, rel_floor)
    cell = SPD.PeriodicSquareUnitCell(n, n)
    psi = SPD.product_pxp_pepskit_state(cell; state = :down, D = 1)
    params = SPD.TrotterParams(dt, 1, :real, maxdim, 0.0; schedule = :serial, rel_floor = rel_floor)
    total = dt * nsteps
    rep = SPD.validate_pxp_reversibility(psi, total; params)
    println("n=$n dt=$dt nsteps=$nsteps maxdim=$maxdim rel_floor=$rel_floor :",
            " density_drift=", rep.density_drift,
            " blockade_drift=", rep.blockade_drift,
            " energy_drift=", rep.energy_drift,
            " fwd_truncerr=", rep.forward_log.max_truncerr,
            " fwd_nsteps=", rep.forward_log.nsteps)
    flush(stdout)
end

# Only D<=2 is feasible: measure_simple's dense star-energy patch has 4^12
# external legs at D=4 (infeasible). Stay at maxdim in {1,2}.
probe(; n = 4, dt = 0.01, nsteps = 1, maxdim = 1, rel_floor = 0.0)
probe(; n = 4, dt = 0.01, nsteps = 2, maxdim = 2, rel_floor = 0.0)
probe(; n = 4, dt = 0.005, nsteps = 2, maxdim = 2, rel_floor = 0.0)
probe(; n = 4, dt = 0.01, nsteps = 3, maxdim = 2, rel_floor = 0.0)
probe(; n = 4, dt = 0.02, nsteps = 2, maxdim = 2, rel_floor = 0.0)
