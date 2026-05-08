# Notes

This folder collects literature and implementation notes for the 2D triangular-lattice PXP ScarFinder PEPS work.

The PEPS code in this repository is internal tooling for ScarFinder, not a standalone tensor-network package. These notes follow that boundary: they focus on constrained triangular PXP evolution, fixed-bond-dimension PEPS projection, blockade diagnostics, and candidate ranking.

## Files

- `literature/`: downloaded PDFs from arXiv and open paper mirrors.
- `extracted/`: local `pdftotext` output used while preparing the notes.
- `literature_review.md`: systematic review of the ScarFinder, PXP/Rydberg, and PEPS/iPEPS update literature.
- `implementation_roadmap.md`: concrete algorithm decisions for this repository.

## Highest-Priority References

1. Ren, Hallam, Ying, Papic, ScarFinder, PRX Quantum 6, 040332, 2025. arXiv:2504.12383. Local PDF: `literature/scarfinder_2504.12383.pdf`.
2. Lin, Calvera, Hsieh, 2D Rydberg PXP scars, PRB 101, 220304(R), 2020. arXiv:2003.04516. Local PDF: `literature/2d_rydberg_scars_2003.04516.pdf`.
3. Jiang, Weng, Xiang, Simple Update, PRL 101, 090603, 2008. arXiv:0806.3719. Local PDF: `literature/simple_update_jiang_0806.3719.pdf`.
4. Dziarmaga, Neighborhood Tensor Update, PRB 104, 094411, 2021. arXiv:2107.06635. Local PDF: `literature/ntu_ipeps_2107.06635.pdf`.
5. Jordan, Orus, Vidal, Verstraete, Cirac, iPEPS, PRL 101, 250602, 2008. arXiv:cond-mat/0703788. Local PDF: `literature/ipeps_jordan_orus_vidal_0703788.pdf`.

## Immediate Takeaway

Implement the algorithm as:

```text
sample PEPS seed
  -> real-time projected triangular PXP evolution
  -> fixed-D PEPS projection by Simple Update
  -> local energy and blockade diagnostics
  -> optional imaginary-time energy correction
  -> candidate ranking
```

Use dense 7-site star gates as the correctness oracle first. Keep the triangular blockade projector explicit and local. Upgrade from Simple Update to NTU only after the dense-star and Simple Update diagnostics are trustworthy.
