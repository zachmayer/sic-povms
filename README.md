# Numerical SIC-POVM fiducials in 108 dimensions, up to d = 39604

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for i ≠ j. Here each SIC is the Weyl–Heisenberg orbit of one fiducial ψ under X^a Z^b. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that a SIC exists in every dimension.

This repo gives numerical fiducials for every d from 197 to 300, and for d = 2208, 5779, 15128 and 39604. Together with the results through 196 of [Bengtsson and Grassl](https://arxiv.org/abs/2512.13201), this covers every d ≤ 300. A search of the cited literature found prior SICs in nine of the 104 dimensions in 197..300:

- In papers: 199, 228 and 259 ([Scott](https://arxiv.org/abs/1703.03993)), and 292 ([Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872)).
- In talks or surveys: 201 ([Grassl, 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 21), and 204, 224, 255 and 288 ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting Grassl).

We found no prior published SIC in the other 95 dimensions. Unpublished work and unindexed data remain outside the search; orbit equivalence with earlier SICs has not been checked. Of the larger dimensions, SICs were known before in 2208, 5779 and 39604 (we found no posted 39604 fiducial); we found none in 15128 ([literature status](RESULTS.md#literature-status)).

There are 109 files in 108 dimensions: 108 standard-Zauner fiducials in `fiducials/` and one Fibonacci antiunitary fiducial in `fibonacci/`. These are numerical solutions, written to 150 digits in [Scott–Grassl's format](https://arxiv.org/abs/0910.5784), with unit norm and ψ_0 real and positive; a real or imaginary part below 1e-180 is written as 0. Each standard file is an eigenvector of Scott–Grassl's U_Z; the construction picks a representative by Galois data, not eigenvalue, so some lie in its ω or ω² eigenspace (ω = e^(2πi/3)).

This repo supplies explicit fiducials. Radchenko and Wheeler's [*Real quadratic fields and finite quantum dilogarithms I*](https://arxiv.org/abs/2609.21892) and [*Stark units for real quadratic fields and reciprocity laws*](https://www.danrad.net/papers/stark-reciprocity.pdf) claim the results needed for Zauner's conjecture in general; we have not checked them.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198, ran the sweep to 300, and wrote `stark.jl` (from Argon's), `search.py`, `refine.py` and this README; GPT-6 Astra found both d = 2208 SICs; GPT-6 wrote the FFT verifier and the large-d pipeline, found d = 5779 and 15128, and fixed the conjugation precision that d = 39604 needed.*

## Checking the files

`verify.py` shares no code with the construction scripts. It checks all d² overlaps, norm, gauge and the stated symmetry, using O(d) working memory. Its header gives the equations and conventions, including the Fibonacci symmetry of [Grassl–Scott](https://arxiv.org/abs/1707.02944v2), Eqs. 24–25.

```sh
uv run verify.py             # all published files; NumPy, threshold 1e-12
uv run verify.py --precise   # FLINT 640-bit balls, threshold 1e-147
```

The precise mode bounds errors in the saved decimal coordinates; these remain numerical checks. Python 3.12+ and uv suffice for the verifier; the large-d numerical Python scripts pin Python 3.13 and their dependencies in their headers. CI runs the default check on every PR and push, and the precise check weekly; 39604 gets the default check weekly and its precise check by hand (about 12 hours). [RESULTS.md](RESULTS.md) records verification measurements and construction history.

## Methods

**d = 198:** `search.py` implements Scott's numerical search in a Zauner eigenspace. `modal_search.py` distributes independent starts across [Modal](https://modal.com) containers; 800 cores found this SIC within an hour. Search times and the first SIC depend on floating-point details.

**The other 103 in 197..300:** `stark.jl` implements the algebraic method of [Appleby, Flammia and Kopp](https://arxiv.org/abs/2501.03970), using Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit. Its comments explain the replacements for costly library steps and its canonical representative rule. `refine.py` takes the resulting seed to the published precision.

**Large dimensions:** `stark_large.jl` builds the companion-form cocycle ghost directly at full precision from a double-sine series, then streams units and traces. `reconstruct_large.py` completes a partial projector with FFT rows; `refine_large.py` refines it in full space. Their headers explain the departures needed to avoid the smaller scripts' dense matrices. Pin rationale is in `representatives.toml`; the Julia header states the unpinned rule tested by `test_large.py`. The 39604 file comes from the prototype these scripts were ported from; its rebuild with them is deferred ([provenance](RESULTS.md#39604-provenance)).

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

On macOS, `tmutil addexclusion <work dir>` keeps the tensors out of Time Machine backups. Local snapshots still pin each deleted tensor predecessor (about 48 GB at d = 57966) until macOS purges them; if the disk fills, the run stops at a checked write, and `sudo tmutil deletelocalsnapshots /` frees the space before resuming.

Construction holds a kernel lock on its work directory, released when the process exits. Ctrl-C during a threaded stage may not exit promptly; `kill -9` is safe to resume from.

```sh
d=5779
work="${LARGE_D_WORK:-large_d_work/d$d}"
LARGE_D_WORK="$work" julia --project=. --threads=8 --gcthreads=1 --heap-size-hint=4G stark_large.jl "$d" units
LARGE_D_WORK="$work" julia --project=. --threads=8 --gcthreads=1 --heap-size-hint=4G stark_large.jl "$d" tensors
uv run reconstruct_large.py "$work/export" --work "$work/reconstruction" --output "fiducial_${d}_stark.txt"
uv run refine_large.py "fiducial_${d}_stark.txt" --output "fiducials/fiducial_${d}.txt" --workers 8
```

`LARGE_D_WORK=<dir>` selects the work directory (default `large_d_work/d<d>`). Construction uses 2048 bits; if the integer-relation guard requests more, append the bit count after the phase in both commands and use a fresh work directory. Unit rows and tensor slices use up to eight Julia threads in separate processes. Tensor values live in fixed-size file records; each trace axis commits before its predecessor is deleted. Refinement puts its checkpoints, shards and `.jsonl` log under `refine_<input-stem>/` there; `--work <dir>` overrides that refinement directory. Each stage checks its own source and inputs. Python-only changes preserve Julia's construction checkpoints. A Julia source, platform or precision mismatch needs a fresh `LARGE_D_WORK`; for a reconstruction or refinement mismatch, move that subdirectory aside or pass `refine_large.py` a fresh `--work`. Otherwise rerun the same commands to resume.

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

The local test covers small odd/even dimensions, precise verification, recovery from a real numerical checkpoint, byte comparisons, negative controls and the conjugation-precision regression (`test_conjugation.jl`). The Modal runner defaults to the four earlier large files, comparing SHA-256 in isolated copies. Its header explains how to retrieve logs and results.

After rebuilding, require `git diff --exit-code fiducials/ fibonacci/`.

## To do

- **January 2027:** rebuild d = 39604 with the repository's scripts on a suitably sized 8-core machine; [acceptance checks](RESULTS.md#39604-provenance).

Public domain (CC0).
