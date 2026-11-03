# SIC-POVMs in every dimension from 197 to 300, plus 2208, 5779 and 15128

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for i ≠ j. Here each SIC is the Weyl–Heisenberg orbit of one fiducial ψ under X^a Z^b. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that a SIC exists in every dimension.

This repo gives numerical fiducials for every d from 197 to 300. Together with the results through 196 of [Bengtsson and Grassl](https://arxiv.org/abs/2512.13201), this covers every d ≤ 300. A search of the cited literature found prior SICs in nine of our 104 dimensions:

- In papers: 199, 228 and 259 ([Scott](https://arxiv.org/abs/1703.03993)), and 292 ([Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872)).
- In talks or surveys: 201 ([Grassl, 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 21), and 204, 224, 255 and 288 ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting Grassl).

We found no prior published SIC in the other 95 dimensions. Unpublished work and unindexed data remain outside the search; orbit equivalence with earlier SICs has not been checked. Of the larger dimensions, SICs were known before in 2208 and 5779; we found none in 15128 ([literature status](RESULTS.md#literature-status)).

There are 108 files in 107 dimensions: 107 standard-Zauner fiducials in `fiducials/` and one Fibonacci antiunitary fiducial in `fibonacci/`. These are numerical solutions, written to 150 digits in [Scott–Grassl's format](https://arxiv.org/abs/0910.5784), with unit norm and ψ_0 real and positive; a real or imaginary part below 1e-180 is written as 0. Each standard file is an eigenvector of Scott–Grassl's U_Z; the construction picks a representative by Galois data, not eigenvalue, so some lie in its ω or ω² eigenspace.

This repo supplies explicit fiducials. [*Stark units for real quadratic fields and reciprocity laws*](https://www.danrad.net/papers/stark-reciprocity.pdf) claims the results needed for Zauner's conjecture in general; we have not checked it.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198, ran the sweep to 300, and wrote `stark.jl` (from Argon's), `search.py`, `refine.py` and this README; GPT-6 Astra found both d = 2208 SICs; GPT-6 wrote the FFT verifier and the large-d pipeline, found d = 5779 and 15128, and fixed the conjugation precision that d = 39604 needed.*

## Checking the files

`verify.py` shares no code with the construction scripts. It checks all d² overlaps, norm, gauge and the stated symmetry, using O(d) working memory. Its header gives the equations and conventions, including the Fibonacci symmetry of [Grassl–Scott](https://arxiv.org/abs/1707.02944v2), Eqs. 24–25.

```sh
uv run verify.py             # all published files; NumPy, threshold 1e-12
uv run verify.py --precise   # FLINT 640-bit balls, threshold 1e-147
```

The precise mode bounds errors in the saved decimal coordinates; these remain numerical checks. Python 3.12+ and uv suffice for the verifier; the large-d numerical Python scripts pin Python 3.13 and their dependencies in their headers. CI runs the default check on every PR and push, and the precise check weekly. [RESULTS.md](RESULTS.md) records verification measurements and construction history.

## Methods

**d = 198:** `search.py` implements Scott's numerical search in a Zauner eigenspace. `modal_search.py` distributes independent starts across [Modal](https://modal.com) containers; 800 cores found this SIC within an hour. Search times and the first SIC depend on floating-point details.

**The other 103 in 197..300:** `stark.jl` implements the algebraic method of [Appleby, Flammia and Kopp](https://arxiv.org/abs/2501.03970), using Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit. Its comments explain the replacements for costly library steps and its canonical representative rule. `refine.py` takes the resulting seed to the published precision.

**Large dimensions:** `stark_large.jl` uses a companion-form cocycle and streamed trace calculations. `reconstruct_large.py` completes a partial projector with FFT rows; `refine_large.py` refines it in full space. Their headers explain the departures needed to avoid the smaller scripts' dense matrices. All three published standard-Zauner large-d files retain an unpublished prototype's integer quotient coordinates and determinant shifts in `representatives.toml`; the Julia header states the rule for unpinned dimensions. `test_large.py` tests that rule at small d.

**Fibonacci 2208:** `fibonacci_2208.py` searches the antiunitary fixed space subject to alignment constraints, with its prototype's ambient random seed fixed in the script. It uses the same large-d refiner.

## Reproducing

Install Julia 1.13.2 and uv, then run from this directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=. stark.jl 197
uv run refine.py fiducial_197_stark.txt
uv run refine.py fiducials/fiducial_198.txt           # re-derive from double precision
uv run search.py 19
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60
```

`STARK_CHECKPOINTS=<existing-dir>` saves the smaller construction's expensive steps. Resume with the same command and directory.

For the large-d companion-form construction:

```sh
d=5779
julia --project=. --threads=1 --gcthreads=1 stark_large.jl "$d"
uv run refine_large.py "fiducial_${d}_stark.txt" --output "fiducials/fiducial_${d}.txt" --workers 8
```

`LARGE_D_WORK=<dir>` selects the work directory (default `large_d_work/d<d>`). Construction uses 2048 bits; if the integer-relation guard requests more, supply the second positional argument and a fresh work directory. Refinement puts its checkpoints, shards and `.jsonl` log under `refine_<input-stem>/` there; `--work <dir>` overrides that refinement directory. Checkpoints bind inputs and settings to the code. Use fresh directories after source or platform changes; otherwise rerun the same commands to resume.

For the published Fibonacci file:

```sh
uv run fibonacci_2208.py --fetch-input
uv run refine_large.py fiducial_2208_fibonacci_seed.txt --output fibonacci/fiducial_2208.txt
```

The [48g source supplement](https://arxiv.org/src/1707.02944v2) is fetched separately. The script extracts only `sicfiducial.48.g.150.txt` and checks its SHA-256. This gitignored external input is not covered by the repo's CC0 license. With a local tarball, replace `--fetch-input` by `--source-archive /path/to/src.tar`; with the checked input present, omit both flags.

Local acceptance and fresh Linux rebuilds:

```sh
uv run test_large.py --work /absolute/path/outside/repo
modal run --detach modal_rebuild.py
modal run --detach modal_rebuild.py --files fiducials/fiducial_5779.txt
```

The local test covers small odd/even dimensions, precise verification, recovery from a real numerical checkpoint, byte comparisons and negative controls. The Modal runner rebuilds all four large published files by default, comparing SHA-256 in isolated copies. Its header explains how to retrieve logs and results.

After rebuilding, require `git diff --exit-code fiducials/ fibonacci/`.

Public domain (CC0).
