# sic-povms

- **Merge level: push to `main`.** Claude may push directly or merge its own PRs.
- The README is the paper. Claims must be checkable from this repo; measurements come from runs recorded in a PR. Call the results numerical, never exact, proved or certified.
- Only `refine.py` and `refine_large.py` write published files. Preserve their bytes, format and symmetry conventions. Seeds and work files stay untracked; external inputs need checksums.
- `verify.py` is independent: it imports no repository scripts, and none imports it. Add checks there for new claims. Every file passes both verifier modes.
- Commands belong in README; per-d facts and measurements belong in RESULTS. Keep representative rules in script headers and quotient pins in `representatives.toml`.
- New dimensions need a RESULTS row and updated README counts. Before claiming novelty, check the sources in RESULTS, Scott 2017, Bengtsson–Grassl, Appleby–Flammia–Kopp and recent arXiv.
- Changes to `stark.jl`, `refine.py`, `search.py`, Manifest or Zauner.jl's revision regenerate the 14 algebraic d in 197..253 and standard 2208, plus targeted dimensions. Preserve `stark.jl`'s canonical rules and independence from random state, run order and Hecke's representation.
- Construction or refinement changes rebuild every affected published file from construction. Changed large-d scripts, quotient pins or shared imports require `modal_rebuild.py` on affected files. A saved seed does not suffice.
- Published bytes must agree across platforms and BLAS/thread settings. New or changed files require a second-platform rebuild, matching SHA-256 and both independent verifier modes before merge. Explain any changed file as a new numerical result.
- Even-d displacement labels and symplectic matrices live modulo 2d: D_(p+du) = ±D_p.
