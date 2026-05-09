# Notes

This folder collects literature and implementation notes for the 2D triangular-lattice PXP ScarFinder PEPS work.

The PEPS code in this repository is internal tooling for ScarFinder, not a standalone tensor-network package. The literature kept here is intentionally narrow: tensor-network algorithms that inform fixed-bond-dimension PEPS evolution, projection, truncation diagnostics, and future environment-aware updates. Scar, Rydberg-platform, and broad model-background papers are not part of this curated Notes bibliography.

## Files

- `literature/`: downloaded PDFs from arXiv and open paper mirrors.
- `extracted/`: local `pdftotext` output used while preparing the notes.
- `literature_review.md`: systematic review of PEPS/iPEPS update algorithms relevant to this repository.
- `implementation_roadmap.md`: concrete algorithm decisions for this repository.
- `current_peps_evolution_solution.md`: current implemented PEPS evolution algorithm and backend status.

## Highest-Priority Algorithm References

1. Jiang, Weng, Xiang, Simple Update, PRL 101, 090603, 2008. arXiv:0806.3719. Local PDF: `literature/simple_update_jiang_0806.3719.pdf`.
2. Dziarmaga, Neighborhood Tensor Update, PRB 104, 094411, 2021. arXiv:2107.06635. Local PDF: `literature/ntu_ipeps_2107.06635.pdf`.
3. Phien, Bengua, Tuan, Corboz, Orus, fast Full Update and gauge fixing, PRB 92, 035142, 2015. arXiv:1503.05345. Local PDF: `literature/fast_full_update_1503.05345.pdf`.
4. Czarnik, Dziarmaga, first-principles iPEPS time evolution, PRB 98, 045110, 2018. arXiv:1804.03872. Local PDF: `literature/first_principles_ipeps_1804.03872.pdf`.
5. Orus, Vidal, CTMRG contraction for iPEPS, PRB 80, 094403, 2009. arXiv:0905.3225. Local PDF: `literature/ctm_ipeps_orus_vidal_0905.3225.pdf`.

## Immediate Takeaway

Implement the algorithm as:

```text
sample PEPS seed
  -> real-time projected triangular PXP evolution at dynamics_maxdim
  -> per-gate PEPS projection by Simple Update
  -> optional hard truncation to scar_maxdim after each projection interval
  -> local energy and blockade diagnostics
  -> optional imaginary-time energy correction
  -> candidate ranking
```

Use dense 7-site star gates as the correctness oracle first. Keep the triangular blockade projector explicit and local. The two-tier dimension split improves the evolve-project map but does not remove the need for per-gate iPEPS projection. Upgrade from Simple Update/ring update to NTU only after the dense-star and Simple Update diagnostics are trustworthy. Treat CTMRG, fast Full Update, and variational update papers as future accuracy infrastructure, not as requirements for the first working local algorithm.
