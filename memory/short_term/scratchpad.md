# Scratchpad

Temporary notes can go here during an active session. Keep this file short and
safe to rewrite.

Current scratchpad:

- Priority now: iPEPS+CTM CPU utilization and tensor-operation efficiency.
- First follow-up experiment should be warmed, single-session CTM timing, not
  parallel first-compile probes.
- Candidate command shape:
  `JULIA_NUM_THREADS=42 SQUAREPXP_CTM_BLAS_THREADS=1 SQUAREPXP_CTM_STRIDED_THREADS=42 SQUAREPXP_CTM_STRIDED_THREADED_MUL=true SQUAREPXP_CTM_PEPSKIT_SCHEDULER=dynamic`.
