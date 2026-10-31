# Numerical results

The 89 dimensions in 199..300 that were not already in this repo: `julia --project=. stark.jl <d>` (Julia 1.13.2, one thread), then `refine.py`, then `verify.py --precise` on the Mac.

- Form: the reduced binary quadratic form ⟨a,b,c⟩ of the principal class `stark.jl` chose; q its conductor; h the class number.
- Mac: 12-core Apple silicon, 96 GB, up to 6 jobs at once. Modal: one Linux x86 container per dimension, about 2.2× slower per job; times there are of the last, checkpoint-resumed attempt, so they leave out the Gauss–Newton steps of earlier preempted attempts (273: about 22 hours from first start to finish).
- Times are wall-clock; peak memory is from `time`. 294 ran its Gauss–Newton on Modal and was finished on the Mac from those checkpoints; its row is the Mac run.
- The early dimensions ran on earlier versions of `stark.jl` (each version regenerated the 14 algebraic files of 197..253 byte for byte). Rerunning 206 and 214 with the final version changed their double-precision output in the 17th digit and left the refined files byte-identical; the other early dimensions were not rerun.
- The last column is the literature search in the README; blank means none found.

| d | form | q | h | U_Z eigenvalue | ran on | stark.jl (s) | peak RSS (GiB) | SIC published before |
|---|---|---|---|---|---|---|---|---|
| 199 | ⟨2,-4,1⟩ | 1 | 1 | 1 | Mac | 589 | 3.9 | paper |
| 200 | ⟨1,-199,1⟩ | 1 | 10 | 1 | Mac | 1866 | 6.2 |  |
| 201 | ⟨1,-200,1⟩ | 3 | 20 | 1 | Mac | 2930 | 5.4 | talk or survey |
| 202 | ⟨1,-201,1⟩ | 1 | 8 | 1 | Mac | 790 | 6.7 |  |
| 204 | ⟨1,-203,1⟩ | 1 | 8 | 1 | Mac | 1383 | 5.6 | talk or survey |
| 205 | ⟨1,-204,1⟩ | 1 | 14 | 1 | Mac | 1203 | 6.5 |  |
| 206 | ⟨23,-69,1⟩ | 1 | 8 | e^(2πi·2/3) | Mac | 1087 | 6.6 |  |
| 208 | ⟨1,-207,1⟩ | 1 | 8 | 1 | Mac | 3387 | 5.4 |  |
| 209 | ⟨1,-208,1⟩ | 1 | 16 | e^(2πi·1/3) | Mac | 1603 | 5.5 |  |
| 210 | ⟨1,-209,1⟩ | 3 | 12 | 1 | Mac | 4920 | 5.9 |  |
| 212 | ⟨1,-211,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 2008 | 7.2 |  |
| 213 | ⟨1,-212,1⟩ | 1 | 24 | 1 | Mac | 2499 | 5.7 |  |
| 214 | ⟨1,-213,1⟩ | 1 | 10 | 1 | Mac | 1672 | 7.8 |  |
| 216 | ⟨1,-215,1⟩ | 1 | 12 | 1 | Mac | 3986 | 5.6 |  |
| 217 | ⟨1,-216,1⟩ | 1 | 18 | 1 | Mac | 4470 | 6.3 |  |
| 218 | ⟨1,-217,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 2272 | 6.6 |  |
| 219 | ⟨55,-110,1⟩ | 3 | 12 | 1 | Mac | 582 | 5.1 |  |
| 220 | ⟨1,-219,1⟩ | 1 | 8 | 1 | Mac | 2397 | 6.8 |  |
| 221 | ⟨1,-220,1⟩ | 1 | 24 | e^(2πi·1/3) | Mac | 3598 | 6.4 |  |
| 222 | ⟨1,-221,1⟩ | 1 | 12 | 1 | Mac | 4012 | 6.1 |  |
| 224 | ⟨1,-15,1⟩ | 1 | 2 | e^(2πi·2/3) | Mac | 257 | 4.3 | talk or survey |
| 225 | ⟨1,-224,1⟩ | 1 | 12 | 1 | Mac | 1571 | 5.9 |  |
| 226 | ⟨1,-225,1⟩ | 1 | 15 | 1 | Mac | 4299 | 8.0 |  |
| 228 | ⟨37,-47,1⟩ | 3 | 3 | 1 | Mac | 365 | 5.1 | paper |
| 229 | ⟨1,-228,1⟩ | 1 | 12 | 1 | Mac | 1275 | 6.4 |  |
| 230 | ⟨1,-229,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 3106 | 6.0 |  |
| 231 | ⟨58,-116,1⟩ | 1 | 16 | 1 | Mac | 1701 | 6.6 |  |
| 232 | ⟨1,-231,1⟩ | 1 | 6 | 1 | Mac | 1654 | 6.4 |  |
| 233 | ⟨26,-78,1⟩ | 1 | 12 | e^(2πi·2/3) | Mac | 457 | 5.2 |  |
| 234 | ⟨1,-233,1⟩ | 1 | 16 | 1 | Modal | 19017 | 4.9 |  |
| 235 | ⟨59,-118,1⟩ | 1 | 6 | 1 | Mac | 380 | 5.0 |  |
| 236 | ⟨1,-235,1⟩ | 1 | 18 | e^(2πi·1/3) | Modal | 15905 | 7.2 |  |
| 237 | ⟨1,-236,1⟩ | 3 | 16 | 1 | Mac | 2782 | 5.4 |  |
| 238 | ⟨1,-237,1⟩ | 1 | 10 | 1 | Mac | 3184 | 6.0 |  |
| 240 | ⟨1,-239,1⟩ | 1 | 16 | 1 | Modal | 8899 | 16.2 |  |
| 242 | ⟨3,-27,1⟩ | 1 | 1 | e^(2πi·2/3) | Mac | 283 | 8.4 |  |
| 244 | ⟨5,-35,1⟩ | 1 | 2 | 1 | Mac | 320 | 5.3 |  |
| 246 | ⟨1,-245,1⟩ | 9 | 18 | 1 | Modal | 8961 | 6.1 |  |
| 247 | ⟨62,-124,1⟩ | 1 | 10 | 1 | Mac | 873 | 5.5 |  |
| 248 | ⟨31,-37,1⟩ | 1 | 2 | 1 | Mac | 431 | 5.5 |  |
| 250 | ⟨1,-249,1⟩ | 1 | 10 | 1 | Mac | 2496 | 10.1 |  |
| 252 | ⟨1,-251,1⟩ | 1 | 16 | 1 | Modal | 41269 | 5.5 |  |
| 254 | ⟨1,-253,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 6687 | 7.1 |  |
| 255 | ⟨1,-16,1⟩ | 3 | 2 | 1 | Mac | 331 | 4.8 | talk or survey |
| 256 | ⟨1,-255,1⟩ | 1 | 12 | 1 | Modal | 3821 | 7.6 |  |
| 257 | ⟨1,-256,1⟩ | 1 | 16 | 1 | Mac | 5797 | 19.6 |  |
| 258 | ⟨1,-257,1⟩ | 1 | 16 | 1 | Modal | 23611 | 6.1 |  |
| 259 | ⟨4,-9,1⟩ | 1 | 2 | 1 | Mac | 646 | 4.8 | paper |
| 260 | ⟨29,-87,1⟩ | 1 | 6 | e^(2πi·2/3) | Mac | 2637 | 5.4 |  |
| 261 | ⟨1,-260,1⟩ | 1 | 40 | 1 | Modal | 60281 | 8.0 |  |
| 262 | ⟨1,-261,1⟩ | 1 | 12 | 1 | Mac | 7310 | 6.7 |  |
| 263 | ⟨66,-132,1⟩ | 1 | 8 | e^(2πi·1/3) | Mac | 480 | 5.6 |  |
| 264 | ⟨1,-263,1⟩ | 3 | 16 | 1 | Modal | 27481 | 8.7 |  |
| 265 | ⟨1,-264,1⟩ | 1 | 16 | 1 | Mac | 2404 | 6.6 |  |
| 266 | ⟨1,-265,1⟩ | 1 | 18 | 1 | Modal | 32328 | 6.6 |  |
| 267 | ⟨67,-134,1⟩ | 1 | 8 | 1 | Mac | 747 | 5.3 |  |
| 268 | ⟨1,-267,1⟩ | 1 | 16 | 1 | Modal | 44851 | 9.4 |  |
| 269 | ⟨30,-90,1⟩ | 1 | 8 | e^(2πi·2/3) | Mac | 404 | 5.7 |  |
| 270 | ⟨1,-269,1⟩ | 1 | 10 | 1 | Mac | 6476 | 7.4 |  |
| 271 | ⟨17,-68,1⟩ | 1 | 4 | 1 | Mac | 375 | 5.5 |  |
| 272 | ⟨1,-271,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 4310 | 7.0 |  |
| 273 | ⟨1,-272,1⟩ | 3 | 36 | 1 | Modal | 38268 | 10.2 |  |
| 274 | ⟨11,-55,1⟩ | 1 | 3 | 1 | Mac | 405 | 11.2 |  |
| 275 | ⟨13,-35,1⟩ | 1 | 2 | e^(2πi·1/3) | Mac | 368 | 5.3 |  |
| 276 | ⟨1,-275,1⟩ | 1 | 16 | 1 | Modal | 35129 | 9.0 |  |
| 277 | ⟨1,-276,1⟩ | 1 | 12 | 1 | Mac | 2713 | 7.6 |  |
| 278 | ⟨5,-19,1⟩ | 1 | 1 | e^(2πi·2/3) | Mac | 455 | 10.4 |  |
| 279 | ⟨70,-140,1⟩ | 1 | 16 | 1 | Mac | 3529 | 7.5 |  |
| 280 | ⟨1,-279,1⟩ | 1 | 10 | 1 | Modal | 15199 | 8.0 |  |
| 281 | ⟨1,-280,1⟩ | 1 | 28 | e^(2πi·1/3) | Mac | 8274 | 8.2 |  |
| 282 | ⟨1,-281,1⟩ | 3 | 14 | 1 | Mac | 9607 | 7.1 |  |
| 283 | ⟨71,-142,1⟩ | 1 | 12 | 1 | Mac | 1358 | 6.8 |  |
| 284 | ⟨1,-283,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 9903 | 8.3 |  |
| 285 | ⟨1,-284,1⟩ | 1 | 16 | 1 | Mac | 6388 | 7.3 |  |
| 286 | ⟨1,-285,1⟩ | 1 | 14 | 1 | Modal | 28612 | 9.5 |  |
| 287 | ⟨2,-24,1⟩ | 1 | 3 | e^(2πi·2/3) | Mac | 353 | 5.5 |  |
| 288 | ⟨1,-17,1⟩ | 1 | 2 | 1 | Mac | 428 | 5.4 | talk or survey |
| 289 | ⟨1,-288,1⟩ | 1 | 24 | 1 | Mac | 7238 | 7.9 |  |
| 290 | ⟨1,-289,1⟩ | 1 | 12 | e^(2πi·1/3) | Mac | 10722 | 7.7 |  |
| 291 | ⟨55,-74,1⟩ | 3 | 8 | 1 | Mac | 528 | 5.7 |  |
| 292 | ⟨17,-19,1⟩ | 1 | 1 | 1 | Mac | 903 | 11.7 | paper |
| 293 | ⟨6,-42,1⟩ | 1 | 4 | 1 | Mac | 454 | 6.2 |  |
| 294 | ⟨1,-293,1⟩ | 1 | 20 | 1 | Modal + Mac | 6775 | 24.2 |  |
| 295 | ⟨74,-148,1⟩ | 1 | 8 | 1 | Mac | 1226 | 6.6 |  |
| 296 | ⟨33,-99,1⟩ | 1 | 10 | e^(2πi·2/3) | Mac | 2253 | 7.6 |  |
| 297 | ⟨37,-44,1⟩ | 1 | 2 | 1 | Mac | 539 | 5.6 |  |
| 298 | ⟨1,-297,1⟩ | 1 | 12 | 1 | Mac | 7921 | 8.0 |  |
| 299 | ⟨3,-30,1⟩ | 1 | 2 | 1 | Mac | 373 | 5.8 |  |
| 300 | ⟨1,-299,1⟩ | 3 | 12 | 1 | Modal | 24952 | 9.1 |  |

## Large-d results

The standard-Zauner files use `stark_large.jl <d>` (2048-bit construction), then `refine_large.py` (640-bit refinement). The additional Fibonacci file uses `fibonacci_2208.py --fetch-input` and the same refiner. All four files are numerical results. Every published large-d file uses data pinned in `representatives.toml` from an unpublished prototype run so its Mac runs and Linux rebuilds compare byte for byte. The pins select integer quotient coordinates, determinant shifts and the Fibonacci search seed; only `test_large.py` at small d tests the canonical rule derived from d.

*This large-d results section is mostly AI-generated by GPT-6.*

| d | standard-Zauner form | conductor q | class number h | U_Z eigenvalue |
|---|---|---:|---:|---|
| 2208 | ⟨1,−7,1⟩ | 3 | 1 | 1 |
| 5779 | ⟨1,−3,1⟩ | 1 | 1 | 1 |
| 15128 | ⟨1,−3,1⟩ | 1 | 1 | exp(4πi/3) |

The separate `fibonacci/fiducial_2208.txt` has F_f symmetry, whose 16th power is type F_a; the standard 2208 file has type F_z ([Grassl–Scott, Proposition III.5](https://arxiv.org/abs/1707.02944v2)).

### Linux rebuild measurements

`modal_rebuild.py` launched each file's construction and refinement commands in a fresh isolated copy on Linux x86, with 8 vCPU, 64 GiB reserved memory and Julia 1.13.2. Julia used one language thread and eight BLAS threads; the Python scripts set one BLAS/FFT thread per process, with eight refinement workers. Python dependencies are pinned in the script headers. The times below are wall-clock seconds recorded by the runner, rounded to two decimal places, and exclude image building. They measure the repository's commands, including reconstruction in the standard construction stage.

| published file | construction (s) | refinement (s) | Linux rebuild |
|---|---:|---:|---|
| `fiducials/fiducial_2208.txt` | 1,288.82 | 286.65 | Byte-identical |
| `fibonacci/fiducial_2208.txt` | 91.62 | 291.13 | Byte-identical |
| `fiducials/fiducial_5779.txt` | 9,634.27 | 4,915.08 | Byte-identical |
| `fiducials/fiducial_15128.txt` | TBD-15128 | TBD-15128 | Running; no completed hash comparison |

The 64 GiB reservation is not a peak-RSS measurement. Peak RSS is unavailable in the retained local records; `/usr/bin/time -v` records it in the stage logs on the `sic-large-rebuild` Volume.

The runner compares SHA-256 of the rebuilt and published bytes. The completed result records report matching hashes:

| published file | run ID | matching SHA-256 |
|---|---|---|
| `fiducials/fiducial_2208.txt` | `2d79f5ee3be840328548c6ea125e33da` | `8d922ad71a48d1c54c887d094bfa8e61f255455ddd1da4e3d4c4ff270f225d48` |
| `fibonacci/fiducial_2208.txt` | `345d62c286d949c595eb1f029d92adc9` | `811ca49245b4c36b0027a0ed42895073e32108104fd15880aeddf57f191eb435` |
| `fiducials/fiducial_5779.txt` | `feb416d67f27421bb3fcadf9bfd01a7b` | `3f5ccd8751146214272e057ed47c97201e40a3e8f481d45ad0f4ffb6cd2ac48a` |

The 15128 run ID is `83d8ca5483764e0b855c7cb21047df57`. Its local result record reports construction still running as of 31 October 2026; its timing cells remain `TBD-15128`. These runs used eight Julia BLAS threads; they do not establish equality across all thread counts or BLAS builds.

### Saved-file verification

`uv run verify.py` checks all 108 published files in about 13 s on an Apple M2 Max, including uv startup. `uv run verify.py --precise` checks the saved decimal coordinates with 640-bit balls; its maximum deviation includes the unit-norm error, and its symmetry residual uses U_Z or F_f according to the folder. Both upper bounds must be below 1e-147 for PASS. These are numerical checks of the coordinates.

A fresh offline default sweep passed all 108 files in 11.96 s, including uv startup, using Python 3.13 with the pinned NumPy 2.5.4 and python-flint 0.9.0 dependencies.

| published file | maximum deviation, precise upper bound | U_Z residual, precise upper bound | U_Z eigenvalue |
|---|---:|---:|---|
| `fiducials/fiducial_5779.txt` | 4.4e-153 | 9.1e-152 | ω⁰ |
| `fiducials/fiducial_15128.txt` | 3.0e-152 | 6.6e-152 | ω² |

Both ran in one `uv run verify.py fiducials/fiducial_5779.txt fiducials/fiducial_15128.txt --precise` call on the M2 Max: 4,504 s wall time in total, while other jobs shared the machine.

### Literature status

- **2208:** Horodecki, Rudnicki and Życzkowski's [2020 survey, §II.1](https://arxiv.org/abs/2002.03233) credits Grassl with a numerical SIC in this dimension. Grassl's [2021 slides, p. 14](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf) also report it. We found no posted 2208 fiducial in the tables described below. Orbit equivalence to either of this repository's 2208 files has not been checked.
- **5779:** a prior SIC is included in Appleby, Bengtsson, Grassl, Harrison and McConnell's [Stark-unit construction](https://arxiv.org/abs/2112.05552), with [5779 data posted by Grassl](http://s208785153.online.de/SICPOVM/5779/). This repository's result is numerical; orbit equivalence has not been checked.
- **15128:** we found no prior SIC in the sources checked. [Bengtsson–McConnell, arXiv:2606.23535v1, §2](https://arxiv.org/html/2606.23535v1#S2), dated 22 June 2026, explicitly states that at d = 15128 “as yet no SIC has been found.” This dates the literature claim; it does not establish the absence of later or unpublished work.

The literature check used these additional records:

- [Grassl's SIC tables](http://s208785153.online.de/SICPOVM/), the frame target of [sicpovm.markus-grassl.de](http://sicpovm.markus-grassl.de/), are reachable over plain HTTP. The checked tables list 5779 with linked data and have no 2208 or 15128 entry.
- Grassl's [UiO seminar abstract for 17 September 2026](https://www.mn.uio.no/math/english/research/projects/qombine/events/seminars/2026-09-17-grassl.html) reports Stark-unit SICs “in some four- and five-digit dimensions, the largest being 39604”. It does not enumerate those dimensions, so this statement alone does not resolve the status of 15128.
- [Bengtsson–Grassl–McConnell, arXiv:2403.02872v1, §9](https://arxiv.org/html/2403.02872v1#S9) reports partially unpublished higher-dimensional constructions. This is another limit on what a search of posted fiducials can establish.

These sources support a report of what we found in the literature, not a claim of novelty. Unpublished work and unindexed data remain outside the search.
