# SIC-POVMs in every dimension from 197 to 300, and in 2208

A SIC (symmetric informationally complete POVM) in dimension d is a set of d² unit vectors in ℂ^d with |⟨ψ_i|ψ_j⟩|² = 1/(d+1) for every pair i ≠ j. Zauner conjectured in 1999 ([thesis](http://www.gerhardzauner.at/documents/gz-quantendesigns.pdf); [English](https://doi.org/10.1142/S0219749911006776)) that one exists in every dimension. SICs are known in every dimension up to 196 ([Bengtsson and Grassl](https://arxiv.org/abs/2512.13201)) and in some higher ones.

This repo has a numerical Weyl–Heisenberg SIC fiducial for every d from 197 to 300 (104 files), so with Bengtsson and Grassl's, every d ≤ 300 has one. Which of these dimensions are new rests on a literature search (Grassl's tables, Scott 2017, Bengtsson–Grassl, Appleby–Flammia–Kopp, recent arXiv), and our view of unpublished and talk-only results is spotty:

- **No SIC found published before (95 dimensions):** every d in 197..300 except the nine below.
- **SIC known before (9):** in papers, 199, 228 and 259 ([Scott](https://arxiv.org/abs/1703.03993)) and 292 ([Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872)); only in talks or surveys, 201 ([Grassl, 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 21) and 204, 224, 255 and 288 ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting Grassl). Ours may lie in different orbits from those; we have not checked orbit equivalence.

It also has one in d = 2208, a Fibonacci–Lucas dimension ([Grassl and Scott](https://arxiv.org/abs/1707.02944)). Grassl has reported a numerical Fibonacci–Lucas SIC there ([Horodecki, Rudnicki and Życzkowski](https://arxiv.org/abs/2002.03233), crediting his 2020 talk; [2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 12); we found no posted fiducial. In d = 2208 a Fibonacci–Lucas SIC's order-3 symmetry is of type F_a, not Zauner's (Grassl and Scott, Proposition III.5), and no SIC orbit has been seen with both types ([Appleby, Flammia and Kopp](https://arxiv.org/abs/2501.03970), Question 2), so ours is presumably in a different orbit from his. One of his type is in `fibonacci/` (below); we have not compared orbits.

This is about explicit fiducials. A recent manuscript, [*Stark units for real quadratic fields and reciprocity laws*](https://www.danrad.net/papers/stark-reciprocity.pdf), claims the results needed for Zauner's conjecture in general; we have not checked it.

Each SIC is the orbit of a fiducial ψ under the d² displacements D_p = X^a Z^b (up to phase), so ψ determines it. `fiducials/fiducial_<d>.txt` gives ψ to 150 digits in the file format and conventions of Scott and Grassl ([arXiv:0910.5784](https://arxiv.org/abs/0910.5784)) and Scott (arXiv:1703.03993): unit norm, ψ_0 real and positive, and an eigenvector of their order-3 Zauner unitary U_Z. These are numerical solutions, not exact ones. `verify.py` checks all of it and shares no code with the other scripts:

```sh
python3 verify.py            # every file in fiducials/, double precision; max deviation at most about 5e-17 (Python 3.10+, standard library)
uv run verify.py --precise   # 160 digits; PASS means deviation and Zauner residual below 1e-147 (here at most 1e-150); minutes per file
```

Both are direct O(d³) sums: at d = 2208 the first takes about 17 minutes on an M2 Max and `--precise` about 17 hours. GitHub Actions runs the first on every pull request and push to `main`, and weekly.

## How they were found

**d = 198** by Scott's numerical method: random starts in the eigenvalue-1 eigenspace of U_Z, each minimized with L-BFGS (`search.py`). It ran on 800 cloud cores with [Modal](https://modal.com) (`modal_search.py`) and found this SIC within its first hour.

**The other 103** by the algebraic method of Appleby, Flammia and Kopp ([arXiv:2501.03970](https://arxiv.org/abs/2501.03970)), which builds a fiducial from Stark units instead of searching for one. `stark.jl` computes a ghost's overlaps on one Galois orbit to thousands of bits (Gauss–Newton), recognizes integral trace forms of them by LLL in a ring class field H of ℚ(√((d+1)(d−3))), applies an automorphism of H that negates the square root, reconstructs the SIC's overlaps on that orbit, and completes the fiducial by alternating projections. It runs Flammia's [Zauner.jl](https://github.com/sflammia/Zauner.jl), pinned by commit, handles odd and even d and forms beyond Zauner.jl's table, and replaces the steps that fail or take hours at these sizes (each replacement says why). One thread per dimension: from 4 minutes to about a day at class number h ≈ 40, under 25 GB. Its choices are canonical (see its header): of a dimension's Galois-conjugate SICs it returns a fixed one, and of that SIC's images under the Clifford unitaries that commute with U_Z, a fixed fiducial. The rule goes by Galois data, not eigenvalue, so some files lie in U_Z's ω or ω² eigenspace (`verify.py` prints each file's).

`refine.py` takes each fiducial from double precision to 150 digits by Gauss–Newton in 640-bit arithmetic, and writes a real or imaginary part below 1e-180 as 0 (the converged residual is about 1e-194).

**d = 2208** by the same method at order conductor q = 3 (form ⟨1,−7,1⟩, discriminant 45, class number 1), in `stark_2208.jl`, which reuses `stark.jl`. It evaluates the Shintani–Faddeev cocycle as a product along the orbit of the form's companion matrix instead of through long q-Pochhammer recurrences. It then streams ghost units and searches determinant classes of shifts. A 640-bit Gauss–Newton refinement like `refine.py`'s, using FFTs (`refine_2208.py`), took it to 150 digits. Since 2208 ≡ 3 mod 9, the conductor fixes the type of the order-3 symmetry ([Appleby, Flammia and Kopp](https://arxiv.org/abs/2501.03970), Theorem 7.21): q = 3 gives Zauner's, and q = 1 gives F_a, the type of Grassl and Scott's Fibonacci–Lucas SIC here. `fibonacci/fiducial_2208.txt` is of type F_a; it came from a numerical search constrained to align with their d = 48 fiducial (48 → 48·46), now in `fibonacci_2208.py` followed by `refine_2208.py`. `python3 verify.py fibonacci/fiducial_2208.txt` finds its overlaps right (max deviation 9.7e-18) and fails it for not being an eigenvector of U_Z; nothing here checks its F_a symmetry.

[RESULTS.md](RESULTS.md) gives, for the 89 dimensions in 199..300 added last, the form, conductor, class number, eigenvalue, machine, run time and memory.

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

For `fiducials/fiducial_2208.txt`:

```sh
julia --project=. --threads=1 stark_2208.jl
uv run refine_2208.py fiducial_2208_stark.txt --output fiducials/fiducial_2208.txt
uv run --no-project python3 verify.py fiducials/fiducial_2208.txt
```

`stark_2208.jl` pins the original quotient generators and reconstruction choices. Its BLAS thread count is fixed at four; run this stage alone when limiting work to four cores.

For `fibonacci/fiducial_2208.txt`:

```sh
uv run fibonacci_2208.py
uv run refine_2208.py fiducial_2208_fibonacci_seed.txt --output fibonacci/fiducial_2208.txt
```

`fibonacci_2208.py` fixes seed 2348 and trial 3, recomputing the symmetry basis and alignment constraints. Its input is the bundled `fiducial_48g.txt`, Grassl–Scott's [48g source supplement](https://arxiv.org/src/1707.02944), SHA-256 `f0b3517ccf6e9ebe40eea6cee62b436f80226e5bf99280a2f0b71f894dc9d090`; the script checks it before use. Both Python scripts pin their package versions and use one numerical thread.

Fresh-copy rebuilds on this Mac, with installed package caches, Julia 1.13.3 and Python 3.13.16:

| File | Construction (s) | Refinement (s) | Peak RSS (GiB) |
|---|---:|---:|---:|
| `fiducials/fiducial_2208.txt` | 572.1 | 1029.4 | 16.00 |
| `fibonacci/fiducial_2208.txt` | 13.5 | 902.3 | 0.153 |

The Julia stage ran alone; both single-thread refinements ran together. Peak RSS is the larger stage peak for each file. Both rebuilds are byte-identical to the committed files:

```sh
shasum -a 256 fiducials/fiducial_2208.txt fibonacci/fiducial_2208.txt
# 26ac55517e5ffe5c3bd290654c0e6a3a682b4f0954e52c33f8a210ec1916b292  fiducials/fiducial_2208.txt
# 8242c58e3fb788d0deab69e10025ebe96a48b831488fa2d99d9072e6aa51235d  fibonacci/fiducial_2208.txt
```

The unchanged default verifier passed the rebuilt standard file in 1117.1 s, with peak RSS 22.9 MiB (Python 3.14.8).

`stark.jl` then `refine.py` rewrite the 14 algebraic files from 197 to 253 byte for byte, whatever the random seeds or the order of dimensions (Julia 1.13.2, macOS on Apple silicon). The early ones among the other 89 ran on earlier versions of `stark.jl` (each of which regenerated those 14 byte for byte); of those, only 206 and 214 were rerun with this one, and their refined files came out byte-identical. Of d = 199..300, 16 ran on Modal (Linux x86, same Project and Manifest), 294 on Modal then the Mac, and the rest on the Mac. The double-precision `stark.jl` output differs between the two in the 17th digit, but the refined files agreed byte for byte in both dimensions run on both (197 and 224). Searches are not repeatable: each trial's path depends on floating-point details, so a repeat finds SICs at its own pace.

*Mostly AI-generated: Gemini 4 Argon first found the algebraic SICs; Claude Opus 5.5 found d = 198, ran the sweep to 300, and wrote the code (`stark.jl` from Argon's) and this README; GPT-6 Astra found both d = 2208 SICs; GPT-6 packaged their reproduction code.*

Public domain (CC0).
