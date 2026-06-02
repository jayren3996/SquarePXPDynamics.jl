# The D>1 simple-density gap is inherent, not a bug

`local_density_simple` (`Observables.jl`) computes the canonical Vidal
single-site mean-field RDM *exactly* (a hand-built λ² RDM matches it to machine
zero). Its ~2e-4 disagreement with ED for D>1 is the inherent single-site
Bethe/mean-field failure on a loopy lattice — NOT a measurement bug:
1-site off ~3.5e-4, but 2-site matches to ~4e-27, 5-site star to ~3e-7.

**Do not "fix" `Observables.jl` for this.** D>1 physics must route through the
trusted observable (`exact_density_finite` / CTM), never the simple/local path.
The `density_error_simple > 1e-4` lower-bound tests correctly DOCUMENT this offset
— keep them; they are not traps.

NOTE (regime caveat): at ultra-short time the offset can vanish if truncation
makes the state effectively a product (e.g. one step at D=2 under rel_floor=1e-3 →
bond dim 1 → simple=exact, offset ~7e-7). The offset-documenting tests therefore
run at an entangled time (t=0.1, offset ~6e-3); see
`../stage2-truncation/rel-floor.md`.

Related: `../methodology/revival-validation.md` (the D-convergence rule — D=1 is
never validation; the simple observable is never the convergence metric).
