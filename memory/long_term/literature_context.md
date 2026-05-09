# Literature Context

- Confirmed: Simple Update is the baseline PEPS projection backend because it is cheap, stable, local, and provides useful discarded-weight or residual diagnostics without requiring CTMRG.
  - Source: `Notes/literature_review.md`
  - Source: `Notes/implementation_roadmap.md`

- Confirmed: Neighborhood Tensor Update is the preferred next algorithmic target after the local Simple Update path and diagnostics are trustworthy. It should share the same high-level star-gate interface.
  - Source: `Notes/literature_review.md`
  - Source: `Notes/implementation_roadmap.md`

- Confirmed: CTMRG, fast Full Update, first-principles iPEPS evolution, environment recycling, and variational real-time PEPS are later accuracy infrastructure, not prerequisites for the first ScarFinder-facing workflow.
  - Source: `Notes/literature_review.md`
  - Source: `Notes/README.md`

- Confirmed: PESS/simplex-aware triangular ansatz literature is relevant background for frustrated triangular lattices, but the repo currently keeps a PEPS convention with one physical leg and six virtual legs per site.
  - Source: `Notes/literature_review.md`

- Confirmed: The current curated Notes bibliography intentionally focuses on tensor-network algorithms relevant to fixed-bond-dimension PEPS evolution and excludes broad scar/Rydberg platform background papers.
  - Source: `Notes/README.md`
