# SIC-POVMs in every dimension from 197 to 300, plus 2208, 5779 and 15128

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones.

This repo has a numerical Weyl–Heisenberg SIC fiducial for every d from 197 to 300 (104 files), so with Bengtsson and Grassl's, every d ≤ 300 has one. Which of these dimensions are new rests on a literature search (Grassl's tables, Scott 2017, Bengtsson–Grassl, Appleby–Flammia–Kopp, recent arXiv), and our view of unpublished and talk-only results is spotty:

- **No SIC found published before (95 dimensions):** every d in 197..300 except the nine below.
- **SIC known before (9):** in papers, 199, 228 and 259 ([Scott](https://arxiv.org/abs/1703.03993)) and 292 ([Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872)); only in talks or surveys, 201 ([Grassl, 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 21) and 204, 224, 255 and 288 ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting Grassl). Ours may lie in different orbits from those; we have not checked orbit equivalence.

There are 108 numerical fiducials in 107 dimensions: 107 standard-Zauner files in `fiducials/` and one Fibonacci antiunitary file in `fibonacci/`. [RESULTS.md](RESULTS.md#large-d-results) records the large-d results and their literature status.

This is about explicit fiducials. A recent manuscript, [*Stark units for real quadratic fields and reciprocity laws*](https://www.danrad.net/papers/stark-reciprocity.pdf), claims the results needed for Zauner's conjecture in general; we have not checked it.

Each SIC is the orbit of a fiducial ψ under the d² displacements D_p = X^a Z^b (up to phase), so ψ determines it. `fiducials/fiducial_<d>.txt` gives ψ to 150 digits in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784)) and Scott (arXiv:1703.03993): unit norm, ψ_0 real and positive, and an eigenvector of their order-3 Zauner unitary U_Z. These are numerical solutions, not exact ones. `verify.py` checks these claims and shares no code with the other scripts. For `fibonacci/`, it checks Grassl–Scott's Fibonacci antiunitary F_f instead of U_Z (arXiv:1707.02944v2, Eqs. 24–25); its 16th power is type F_a at d = 2208 (Proposition III.5). Both folders are checked by default:

```sh
uv run verify.py             # every file in fiducials/ and fibonacci/, NumPy double precision; PASS threshold 1e-12
uv run verify.py --precise   # python-flint, 640-bit balls (~192 decimal digits); PASS threshold 1e-147
```

Both modes check all d² overlaps using one FFT row at a time, in O(d² log d) time and O(d) working memory, and apply the symmetry with one Fourier transform and a chirp. Python 3.12+ and uv install the pinned NumPy and python-flint dependencies from the script header. In `--precise`, the maximum deviation (including unit norm) and symmetry residual are computed as rigorous upper bounds for the decimal coordinates in the file; PASS compares the bounds before display rounding. These remain numerical checks of the published coordinates. Measured on an M2 Max, including uv startup: d = 2208 takes 0.27 s in default mode and 107.2 s with `--precise`; the default run over all 106 files takes 0.83 s. GitHub Actions runs the default mode on every file for each pull request and push to `main`, and `--precise` on every file weekly.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198, ran the sweep to 300, and wrote the code (`stark.jl` from Argon's) and this README; GPT-6 Astra found both d = 2208 SICs; GPT-6 packaged their reproduction code and revised its determinism and verification checks; GPT-6 replaced the independent verifier with FFTs; GPT-6 generalized the large-d pipeline, added two numerical results and updated this README.*

## How they were found

**d = 198** by Scott's numerical method: random starts in the eigenvalue-1 eigenspace of U_Z, each minimized with L-BFGS (`search.py`). It ran on 800 cloud cores with [Modal](https://modal.com) (`modal_search.py`) and found this SIC within its first hour.

**The other 103 in 197..300** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one. `stark.jl` computes a ghost's overlaps on one Galois orbit to thousands of bits (Gauss–Newton), recognizes integral trace forms of them by LLL in a ring class field H of ℚ(√((d+1)(d−3))), applies an automorphism of H that negates the square root, reconstructs the SIC's overlaps on that orbit, and completes the fiducial by alternating projections. It runs Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit, handles odd and even d and forms beyond Zauner.jl's table, and replaces the steps that fail or take hours at these sizes (each replacement says why). One thread per dimension: from 4 minutes to about a day at class number h ≈ 40, under 25 GB. Its choices are canonical (see its header): of a dimension's Galois-conjugate SICs it returns a fixed one, and of that SIC's images under the Clifford unitaries that commute with U_Z, a fixed fiducial. The rule goes by Galois data, not eigenvalue, so some files lie in U_Z's ω or ω² eigenspace (`verify.py` prints each file's).

`refine.py` takes each fiducial from double precision to 150 digits by Gauss–Newton in 640-bit arithmetic, and writes a real or imaginary part below 1e-180 as 0 (the converged residual is about 1e-194).

The general large-d pipeline (`stark_large.jl`, `reconstruct_large.py`, `refine_large.py`) factors the Shintani–Faddeev cocycle along the companion-matrix orbit, corrects the ghost with FLINT FFT residuals and a frozen Float64 QR, streams its units and normal-basis traces, then reconstructs from determinant classes without a dense complex d×d matrix. Its header states the canonical choices computed from d; `representatives.toml` preserves the integer quotient coordinates and determinant shifts of the existing representatives. Full-space 640-bit FFT Gauss–Newton refinement writes 150 digits and clips parts below 1e-180 to zero. The separate Fibonacci search uses ambient seeded starts projected onto its antiunitary fixed space and a separately fetched, checksummed Grassl–Scott alignment input; it uses the same general refiner.

[RESULTS.md](RESULTS.md) records forms, conductors, class numbers, eigenvalues, hardware, timings, memory and literature status.

## Reproducing

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once; Julia 1.13.2
julia --project=. stark.jl 197                        # fiducial_197_stark.txt, to double precision
uv run refine.py fiducial_197_stark.txt               # rewrites fiducials/fiducial_197.txt, unchanged
uv run refine.py fiducials/fiducial_198.txt           # re-derives its 150 digits from the first 16, unchanged
uv run search.py 19                                   # Scott's search: a SIC in d = 19 in a second
uvx modal run modal_search.py --d 198 --containers 25 --minutes 60  # a search the size of the one that found d = 198
```

With `STARK_CHECKPOINTS=<dir>` (an existing directory, fresh for each version of `stark.jl`), `stark.jl` saves each Gauss–Newton step and its other expensive intermediates there, and a rerun resumes from them.

For any supported large-d companion form (set `d` to the dimension):

```sh
d=5779
julia --project=. --threads=1 stark_large.jl "$d"
uv run refine_large.py "fiducial_${d}_stark.txt" --output "fiducials/fiducial_${d}.txt" --workers 8
```

Construction defaults to 2048 bits; an optional second positional argument raises it if the integer-relation guard requests more. `LARGE_D_WORK=<dir>` selects its checkpoint directory (default `large_d_work/d<d>`). Refinement checkpoints are bound to the input, code, precision and worker count. Resume with the same commands; use fresh directories after source or platform changes. Each published path has one writer: `refine_large.py` for this pipeline, `refine.py` for the smaller files.

The out7 prototype's measured numerical-stage sums were 1.176–9.776 hours including the independent saved-file check, on an Apple M2 Max with 96 GiB RAM and up to eight cores. Construction accounted for 1,628.79–12,180.58 seconds; refinement for 2,485.15–21,760.65 seconds. The largest recorded single-process Julia peak RSS was 5.66–20.19 decimal GB; multiworker aggregate observations were about 23–32.85 GB and are not continuous peak bounds. These are **out7 prototype timings until the Modal rebuild reports**, not measurements of these packaged commands. The packaged constructor streams unit rows in one Julia process; the prototype used eight unit workers, so that stage will take longer here. Per-file measurements and their limitations are in [RESULTS.md](RESULTS.md#large-d-results).

For the published Fibonacci file:

```sh
uv run fibonacci_2208.py --fetch-input
uv run refine_large.py fiducial_2208_fibonacci_seed.txt --output fibonacci/fiducial_2208.txt
```

Its Grassl–Scott [48g source supplement](https://arxiv.org/src/1707.02944v2) is fetched separately. Only `sicfiducial.48.g.150.txt` is extracted, and its SHA-256 must be `f0b3517ccf6e9ebe40eea6cee62b436f80226e5bf99280a2f0b71f894dc9d090`. The fetched input is gitignored and is not licensed by this repo's CC0. With a local tarball, replace `--fetch-input` by `--source-archive /path/to/src.tar`; with the checked input present, omit both flags.

Fresh Linux rebuilds of selected published paths run independently on Modal:

```sh
modal run --detach modal_rebuild.py
modal run --detach modal_rebuild.py --files fiducials/fiducial_5779.txt
```

The runner uploads this repository snapshot, removes each selected output in an isolated copy, runs its construction and refinement commands, and records logs, stage times, `/usr/bin/time -v` peak RSS and committed/rebuilt SHA-256 on the `sic-large-rebuild` Volume. It prints a run ID and download command. Calls longer than 22 hours continue sequentially from checkpoints, up to eight slices. Each app has one remote function; each selected file has one active writer. Set `SIC_REBUILD_BLAS_THREADS=1` or `8` for separate determinism runs. The default memory reservation is 64 GiB per file; `SIC_REBUILD_MEMORY_MIB=262144` prepares a future larger rebuild. The image build and remote execution still need validation. Local acceptance makes no cloud calls:

```sh
uv run test_large.py --work /absolute/path/outside/repo
```

It constructs small odd/even examples including q = 3, checks the saved digits with the independent precise verifier, and tests the Modal worker's actual local logic, interruption, continuation, hash mismatch and child failure. Before merge, rebuild all affected published files on Linux with matching hashes, then independently verify them at both precisions. After rebuilding in a Git checkout:

```sh
git diff --exit-code fiducials/ fibonacci/
```

`stark.jl` then `refine.py` rewrite the 14 algebraic files from 197 to 253 byte for byte, whatever the random seeds or the order of dimensions (Julia 1.13.2, macOS on Apple silicon). The early ones among the other 89 ran on earlier versions of `stark.jl` (each of which regenerated those 14 byte for byte); of those, only 206 and 214 were rerun with this one, and their refined files came out byte-identical. Of d = 199..300, 16 ran on Modal (Linux x86, same Project and Manifest), 294 on Modal then the Mac, and the rest on the Mac. The double-precision `stark.jl` output differs between the two in the 17th digit, but the refined files agreed byte for byte in both dimensions run on both (197 and 224). The replaced dimension-specific scripts had also passed a Linux rebuild; that evidence does not validate the general scripts. `search.py`'s unconstrained searches depend on floating-point details: a repeat finds SICs at its own pace.

Public domain (CC0).
