# Literature Context

## PXP Scars And ScarFinder

- Confirmed: The project motivation is square-lattice PXP dynamics and
  scar-like candidate search using low-entanglement PEPS/iPEPS methods.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Literature anchors recorded there include ScarFinder
  `https://arxiv.org/abs/2504.12383`, Bernien et al.
  `https://arxiv.org/abs/1707.04344`, Turner et al.
  `https://www.nature.com/articles/s41567-018-0137-5`, and deformed PXP /
  emergent SU(2) work
  `https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.122.220603`.

## iPEPS, CTMRG, And Update Algorithms

- Confirmed: Simple update is the local first algorithmic layer; CTMRG and
  full-update/gauge-fixed truncation are later accuracy infrastructure.
- Source: `notes/2026-05-15-ipeps-literature-code-algorithm-notes.md`
- Literature anchors recorded there include iPEPS gauge fixing / fast full
  update `https://arxiv.org/abs/1503.05345`, CTMRG for iPEPS contraction
  `https://arxiv.org/abs/0905.3225`, and Neighborhood Tensor Update
  `https://arxiv.org/abs/2107.06635`.

## Codebase References

- Confirmed: ITensors.jl is used as the low-level named-index tensor engine.
- Confirmed: PEPSKit.jl is used for experimental CTMRG measurement support,
  but not as the source of the custom five-site PXP update.
- Source: `Project.toml`
- Source: `notes/2026-05-15-pepskit-backend-feasibility.md`
- Source: `src/PEPSKitMeasurements.jl`

## TFIM Benchmark References

- Confirmed: TFIM was chosen as the first non-PXP benchmark target because it
  has published infinite-iPEPS and dynamics references.
- Source: `docs/superpowers/specs/2026-05-15-infinite-tfim-benchmark-design.md`
- Literature anchors in the design include:
  Jordan et al., PRL 101, 250602 (2008),
  `https://doi.org/10.1103/PhysRevLett.101.250602`;
  Schmitt et al., Science Advances 8, eabl6850 (2022),
  `https://arxiv.org/abs/2106.09046`;
  Arias Espinoza and Corboz, PRB 110, 094314 (2024),
  `https://arxiv.org/abs/2405.10628`;
  and Vovrosh et al., finite-system cross-method dynamics context,
  `https://arxiv.org/abs/2511.19340`.
