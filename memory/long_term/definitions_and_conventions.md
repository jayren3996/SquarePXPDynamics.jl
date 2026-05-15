# Definitions And Conventions

## Basis And Operators

- Confirmed: Physical basis convention is `|up> = |0>`, `|down> = |1>`.
- Confirmed: In dense/ITensor code, basis index `1` corresponds to `:up` and
  basis index `2` corresponds to `:down`.
- Confirmed: `:up` is the Rydberg/excited state and `:down` is the
  vacancy/unexcited state in the PXP convention.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`

## Coordinates And Directions

- Confirmed: Square coordinates are `(x, y)`.
- Confirmed: Direction order is `:right`, `:up`, `:left`, `:down`.
- Confirmed: Dense square-star site order is `(center, right, up, left, down)`.
- Source: `README.md`
- Source: `notes/README.md`
- Source: `src/SquareGeometry.jl`
- Source: `src/SquareUnitCells.jl`

## Unit Cells And Scheduling

- Confirmed: Periodic iPEPS work uses rectangular unit cells represented by
  one-based `SquareCoord(x, y)` representatives.
- Confirmed: The deterministic five-color schedule uses
  `color(x, y) = mod(x + 2y, 5) + 1`.
- Confirmed: Five-color compatible periodic cells require dimensions compatible
  with disjoint same-color radius-1 stars; `10 x 10` is the robust default
  referenced in notes and TFIM plans.
- Source: `notes/README.md`
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `src/SquareGeometry.jl`
- Source: `src/SquareUnitCells.jl`

## Tensor And Gate Ordering

- Confirmed: Finite PEPS tensor index order is physical, left, right, up, down.
- Confirmed: Dense five-site star gates are converted to ITensors with primed
  output physical indices and unprimed input physical indices in square-star
  site order.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `src/SquarePEPS.jl`
- Source: `src/SquareIPEPS.jl`

## Public API Discipline

- Confirmed: Exported Julia symbols require docstrings because
  `test/test_public_docs.jl` checks public documentation coverage.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Source: `test/test_public_docs.jl`

## TFIM Convention

- Confirmed: The TFIM benchmark design pins `Z_up_is_plus_one` and an `X`
  transverse field convention for the new star-model layer.
- Confirmed: TFIM star order remains `(center, right, up, left, down)`.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Source: `docs/superpowers/notes/2026-05-15-current-work-infinite-tfim-benchmark.md`
