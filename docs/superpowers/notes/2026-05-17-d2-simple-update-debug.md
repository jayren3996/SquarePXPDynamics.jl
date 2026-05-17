# D2 Simple-Update Debug Note

Date: 2026-05-17

## Question

Does the first M2 audit's D=2 failure come from starting validation with a
padded product state (`maxdim = 2`, link weights like `[1, 0]`) rather than a
minimal product state allowed to grow to update cap D=2?

## Diagnostic

On the same `3 x 3`, all-down, serial PXP setup used by the no-CTM audit
(`total_time = 0.02`, `dt = 0.02`, `cutoff = 1e-12`), compare:

- A: `product_square_ipeps(...; maxdim = 1)`, update cap `maxdim = 1`
- B: `product_square_ipeps(...; maxdim = 1)`, update cap `maxdim = 2`
- C: `product_square_ipeps(...; maxdim = 2)`, update cap `maxdim = 2`

## Result

Case B and Case C agree to roundoff:

| Case | Density | log_norm | max_truncerr |
| --- | ---: | ---: | ---: |
| A | `0.0003989898795957303` | `-1.0773451054821087e-5` | `1.5989338961408006e-7` |
| B | `0.00021094978264193137` | `7.278045395879476` | `1.8551396277413272e-28` |
| C | `0.0002109497826419281` | `7.2780453958794435` | `1.0579672926369623e-28` |

The current validation path therefore does have a padded-product construction,
but that is not the sole cause of the D=2 anomaly. The grow-on-demand D=2 path
already reproduces it.

Per-star tracing shows the first star is fine and divergence begins once
overlapping serial stars grow active D=2 links. For D=2, per-star
normalization increments jump to exact factors such as `sqrt(2)`, `2`, and
`2sqrt(2)`, while truncation remains near zero and reversibility remains small.

An exact dense serial-star circuit for the same `3 x 3` order gives density
`0.00039962698926202146`, close to the D=1 iPEPS value and far from the D=2
simple-observable value. This confirms the issue is not just ED-vs-Trotter
comparison.

## Negative Result

Removing the split scalar from the remaining SVD core reduced the final D=2
`log_norm` from about `7.28` to about `1.04`, but did not change the D=2
density (`0.00021094978264193`). That scalar recycling is therefore a log-norm
symptom, not a sufficient observable fix.

## Decision

Do not spend the next branch on padded-product pruning alone. The next
debugging target is the D=2 QR/SVD star split and/or the D>1 simple-observable
local environment after overlapping serial updates. CTM Stage 2 remains gated
until this pre-CTM D=2 behavior is understood.
