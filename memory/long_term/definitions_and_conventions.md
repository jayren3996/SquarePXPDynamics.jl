# Definitions And Conventions

- Confirmed: Physical basis convention is `|up> = |0>` and `|down> = |1>`.
  - Source: `AGENTS.md`
  - Source: `Notes/implemented_peps_algorithm_detail.md`
  - Source: `src/States.jl`

- Confirmed: Dense 7-site star ordering is center first, followed by the six nearest neighbors in triangular direction order.
  - Source: `AGENTS.md`
  - Source: `Notes/implemented_peps_algorithm_detail.md`
  - Source: `src/Geometry.jl`

- Confirmed: Triangular axial directions are:
  - `1: (1, 0)`
  - `2: (0, 1)`
  - `3: (-1, 1)`
  - `4: (-1, 0)`
  - `5: (0, -1)`
  - `6: (1, -1)`
  - Source: `src/Geometry.jl`

- Confirmed: Opposite virtual directions are `1 <-> 4`, `2 <-> 5`, and `3 <-> 6`.
  - Source: `src/States.jl`

- Confirmed: Dense full-star blockade checks 12 local triangular-star edges: six center-neighbor edges and six neighbor-ring edges.
  - Source: `Notes/implemented_peps_algorithm_detail.md`
  - Source: `src/Models.jl`

- Confirmed: Supported unit cells are `OneSiteUnitCell`, `ThreeSiteUnitCell`, and `SevenSiteUnitCell`; seven-site wrapping follows the star-color function.
  - Source: `src/States.jl`
  - Source: `test/test_states.jl`

- Confirmed: The 7-color schedule uses `star_color(c) = (c.q + 3c.r) mod 7 + 1`, with canonical center `Coord(color - 1, 0)`.
  - Source: `src/Geometry.jl`
  - Source: `src/Evolution.jl`
