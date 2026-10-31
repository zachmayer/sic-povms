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

These numerical files use `stark_large.jl <d>` and `refine_large.py`; the additional Fibonacci file uses `fibonacci_2208.py` and the same refiner. All times below are wall-clock seconds. The two larger rows are **out7 prototype measurements until the packaged Modal rebuild reports**, with eight unit/gradient workers and four BLAS threads for ghost QR. Their stage sums are not uninterrupted stopwatches. Python was 3.13, python-flint 0.9.0, NumPy 2.3.2 and SciPy 1.16.1. The earlier packaging row used NumPy 2.5.4 and SciPy 1.18.1; these are the general Python scripts' pinned versions. Julia construction precision was 2048 bits, refinement 640 bits.

*This large-d results section is mostly AI-generated by GPT-6.*

| d | form | q | h | U_Z eigenvalue | machine | construction / refinement / check (s) | memory | literature status |
|---|---|---:|---:|---|---|---|---|---|
| 2208 | ⟨1,−7,1⟩ | 3 | 1 | 1 | Apple M2 Max, 96 GB; Julia 1.13.3, Python 3.13.16, one numerical thread | Standard: 521.3 / 851.8; Fibonacci search: 61.5 / 908.1 (previous scripts, fresh copies) | Standard peak RSS 15.38 GiB; Fibonacci peak 0.156 GiB | A numerical Fibonacci–Lucas SIC was reported by Grassl ([2020 survey](https://arxiv.org/abs/2002.03233), [2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), p. 12); no posted fiducial found. It has type F_a, whereas this row's standard file has Zauner type F_z ([Grassl–Scott](https://arxiv.org/abs/1707.02944), Prop. III.5; [AFK](https://arxiv.org/abs/2501.03970), Thm. 7.21 and Question 2). The additional `fibonacci/` file has F_f symmetry and type F_a; orbit equivalence to Grassl's was not checked. |
| 5779 | ⟨1,−3,1⟩ | 1 | 1 | 1 | Apple M2 Max, 96 GiB; Julia 1.13.3, up to eight cores | 1,628.79 / 2,485.15 / 119.00; sum 4,232.94 | Ghost peak RSS 5.52 decimal GB; conjugation peak 5.66 GB; unit-worker aggregate about 23 GB (snapshot, not peak bound) | SIC published before: Appleby, Bengtsson, Grassl, Harrison and McConnell, [arXiv:2112.05552](https://arxiv.org/abs/2112.05552), exact literature construction, with data on [Grassl's site](https://sicpovm.markus-grassl.de/). This repository's result is numerical; orbit equivalence was not checked. |
| 15128 | ⟨1,−3,1⟩ | 1 | 1 | exp(4πi/3) | Apple M2 Max, 96 GiB; Julia 1.13.3, up to eight cores | 12,180.58 / 21,760.65 / 1,252.29; sum 35,193.52 | Ghost peak RSS 20.19 decimal GB; conjugation peak 15.44 GB; unit-worker sampled aggregate maximum 32.85 GB (not continuous peak bound) | No published SIC found in the search below; this is a search result, not a novelty proof. |

The last row's driver invocations totaled 9.815 process-timed elapsed hours, over 10.122 UTC hours from first launch to final check. The machine slept on battery during refinement; the difference also includes a gap between invocations. The individual numerical-stage sum excludes some startup and cache-load work. The sum of its eight unit workers' individual RSS maxima was 33.15 decimal GB; those maxima need not coincide and exclude the parent. Neither aggregate observation is a rigorous peak bound.

Out7's independent saved-file scans covered all nonidentity displacements at double precision and 600 bits, before the final zero clipping. Maximum 600-bit squared-overlap errors were 4.36e-153 and 1.90e-153 for the two larger files, respectively; U_Z residuals were 9.09e-152 and 6.65e-152. The first new file is byte-identical to out7. The second differs only at zero-based component 7564, whose imaginary roundoff part (about −3.71e-192) is now zero. The general scripts still need fresh Modal rebuilds and independent saved-file verification before merge.

Literature search for d = 15128 (2026-10-30 session):

- Attempted [Grassl's SIC tables](https://sicpovm.markus-grassl.de/) and the older [SIC URL](https://www.markus-grassl.de/SIC/); the browsing tool could not access the data host. Site-restricted searches for `15128` and `5779` returned no indexed data rows. This does not establish absence from the live tables.
- Read [Bengtsson–Grassl–McConnell, arXiv:2403.02872](https://arxiv.org/abs/2403.02872): `15128` occurs in the dimension-tower sequence, not as a constructed result; checked the explicit dimension lists and Tables 1–3.
- Searched arXiv for `"15128" SIC`, including results since 2024, and for recent `SIC-POVM` work in 2026. The dimension-tower papers [1701.05200](https://arxiv.org/abs/1701.05200) and [1707.09911](https://arxiv.org/abs/1707.09911) discuss the sequence, not an explicit fiducial there. Recent results included [2609.04431](https://arxiv.org/abs/2609.04431), on tower structure. No explicit d = 15128 result was found. Unpublished results and unindexed data remain outside this search.

