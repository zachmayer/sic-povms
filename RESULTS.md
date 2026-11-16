# Numerical results

Construction measurements for 89 of the 104 dimensions in 197..300, each built by `stark.jl` (Julia 1.13.2, one thread) and `refine.py`. The other 15 are d = 198, found by `search.py`, and 14 algebraic files regenerated as described below.

- Form: the reduced binary quadratic form ⟨a,b,c⟩ of the principal class `stark.jl` chose; q its conductor; h the class number; ω = e^(2πi/3).
- Mac: 12-core Apple silicon, 96 GB, up to 6 jobs at once. Modal: one Linux x86 container per dimension, about 2.2× slower per job; times there are of the last, checkpoint-resumed attempt, so they leave out the Gauss–Newton steps of earlier preempted attempts.
- Times are wall-clock; peak memory is from `time`. 294 ran its Gauss–Newton on Modal and was finished on the Mac from those checkpoints; its row is the Mac run.
- The final `stark.jl` and `refine.py` regenerate the 14 algebraic files of 197..253 byte for byte on Julia 1.13.2, macOS Apple silicon, independent of random seeds and dimension order.
- Of the other 89 algebraic files, only 206 and 214 were rerun with the final version; their double-precision seeds changed in the 17th digit and their refined files stayed byte-identical.
- Refined 197 and 224 agreed byte for byte across Mac and Linux, despite double-precision seed differences.
- The last column is the literature search in the README; blank means none found.

| d | form | q | h | U_Z eigenvalue | ran on | stark.jl (s) | peak RSS (GiB) | SIC published before |
|---|---|---|---|---|---|---|---|---|
| 199 | ⟨2,-4,1⟩ | 1 | 1 | 1 | Mac | 589 | 3.9 | paper |
| 200 | ⟨1,-199,1⟩ | 1 | 10 | 1 | Mac | 1866 | 6.2 |  |
| 201 | ⟨1,-200,1⟩ | 3 | 20 | 1 | Mac | 2930 | 5.4 | talk or survey |
| 202 | ⟨1,-201,1⟩ | 1 | 8 | 1 | Mac | 790 | 6.7 |  |
| 204 | ⟨1,-203,1⟩ | 1 | 8 | 1 | Mac | 1383 | 5.6 | talk or survey |
| 205 | ⟨1,-204,1⟩ | 1 | 14 | 1 | Mac | 1203 | 6.5 |  |
| 206 | ⟨23,-69,1⟩ | 1 | 8 | ω² | Mac | 1087 | 6.6 |  |
| 208 | ⟨1,-207,1⟩ | 1 | 8 | 1 | Mac | 3387 | 5.4 |  |
| 209 | ⟨1,-208,1⟩ | 1 | 16 | ω | Mac | 1603 | 5.5 |  |
| 210 | ⟨1,-209,1⟩ | 3 | 12 | 1 | Mac | 4920 | 5.9 |  |
| 212 | ⟨1,-211,1⟩ | 1 | 12 | ω | Mac | 2008 | 7.2 |  |
| 213 | ⟨1,-212,1⟩ | 1 | 24 | 1 | Mac | 2499 | 5.7 |  |
| 214 | ⟨1,-213,1⟩ | 1 | 10 | 1 | Mac | 1672 | 7.8 |  |
| 216 | ⟨1,-215,1⟩ | 1 | 12 | 1 | Mac | 3986 | 5.6 |  |
| 217 | ⟨1,-216,1⟩ | 1 | 18 | 1 | Mac | 4470 | 6.3 |  |
| 218 | ⟨1,-217,1⟩ | 1 | 12 | ω | Mac | 2272 | 6.6 |  |
| 219 | ⟨55,-110,1⟩ | 3 | 12 | 1 | Mac | 582 | 5.1 |  |
| 220 | ⟨1,-219,1⟩ | 1 | 8 | 1 | Mac | 2397 | 6.8 |  |
| 221 | ⟨1,-220,1⟩ | 1 | 24 | ω | Mac | 3598 | 6.4 |  |
| 222 | ⟨1,-221,1⟩ | 1 | 12 | 1 | Mac | 4012 | 6.1 |  |
| 224 | ⟨1,-15,1⟩ | 1 | 2 | ω² | Mac | 257 | 4.3 | talk or survey |
| 225 | ⟨1,-224,1⟩ | 1 | 12 | 1 | Mac | 1571 | 5.9 |  |
| 226 | ⟨1,-225,1⟩ | 1 | 15 | 1 | Mac | 4299 | 8.0 |  |
| 228 | ⟨37,-47,1⟩ | 3 | 3 | 1 | Mac | 365 | 5.1 | paper |
| 229 | ⟨1,-228,1⟩ | 1 | 12 | 1 | Mac | 1275 | 6.4 |  |
| 230 | ⟨1,-229,1⟩ | 1 | 12 | ω | Mac | 3106 | 6.0 |  |
| 231 | ⟨58,-116,1⟩ | 1 | 16 | 1 | Mac | 1701 | 6.6 |  |
| 232 | ⟨1,-231,1⟩ | 1 | 6 | 1 | Mac | 1654 | 6.4 |  |
| 233 | ⟨26,-78,1⟩ | 1 | 12 | ω² | Mac | 457 | 5.2 |  |
| 234 | ⟨1,-233,1⟩ | 1 | 16 | 1 | Modal | 19017 | 4.9 |  |
| 235 | ⟨59,-118,1⟩ | 1 | 6 | 1 | Mac | 380 | 5.0 |  |
| 236 | ⟨1,-235,1⟩ | 1 | 18 | ω | Modal | 15905 | 7.2 |  |
| 237 | ⟨1,-236,1⟩ | 3 | 16 | 1 | Mac | 2782 | 5.4 |  |
| 238 | ⟨1,-237,1⟩ | 1 | 10 | 1 | Mac | 3184 | 6.0 |  |
| 240 | ⟨1,-239,1⟩ | 1 | 16 | 1 | Modal | 8899 | 16.2 |  |
| 242 | ⟨3,-27,1⟩ | 1 | 1 | ω² | Mac | 283 | 8.4 |  |
| 244 | ⟨5,-35,1⟩ | 1 | 2 | 1 | Mac | 320 | 5.3 |  |
| 246 | ⟨1,-245,1⟩ | 9 | 18 | 1 | Modal | 8961 | 6.1 |  |
| 247 | ⟨62,-124,1⟩ | 1 | 10 | 1 | Mac | 873 | 5.5 |  |
| 248 | ⟨31,-37,1⟩ | 1 | 2 | 1 | Mac | 431 | 5.5 |  |
| 250 | ⟨1,-249,1⟩ | 1 | 10 | 1 | Mac | 2496 | 10.1 |  |
| 252 | ⟨1,-251,1⟩ | 1 | 16 | 1 | Modal | 41269 | 5.5 |  |
| 254 | ⟨1,-253,1⟩ | 1 | 12 | ω | Mac | 6687 | 7.1 |  |
| 255 | ⟨1,-16,1⟩ | 3 | 2 | 1 | Mac | 331 | 4.8 | talk or survey |
| 256 | ⟨1,-255,1⟩ | 1 | 12 | 1 | Modal | 3821 | 7.6 |  |
| 257 | ⟨1,-256,1⟩ | 1 | 16 | 1 | Mac | 5797 | 19.6 |  |
| 258 | ⟨1,-257,1⟩ | 1 | 16 | 1 | Modal | 23611 | 6.1 |  |
| 259 | ⟨4,-9,1⟩ | 1 | 2 | 1 | Mac | 646 | 4.8 | paper |
| 260 | ⟨29,-87,1⟩ | 1 | 6 | ω² | Mac | 2637 | 5.4 |  |
| 261 | ⟨1,-260,1⟩ | 1 | 40 | 1 | Modal | 60281 | 8.0 |  |
| 262 | ⟨1,-261,1⟩ | 1 | 12 | 1 | Mac | 7310 | 6.7 |  |
| 263 | ⟨66,-132,1⟩ | 1 | 8 | ω | Mac | 480 | 5.6 |  |
| 264 | ⟨1,-263,1⟩ | 3 | 16 | 1 | Modal | 27481 | 8.7 |  |
| 265 | ⟨1,-264,1⟩ | 1 | 16 | 1 | Mac | 2404 | 6.6 |  |
| 266 | ⟨1,-265,1⟩ | 1 | 18 | 1 | Modal | 32328 | 6.6 |  |
| 267 | ⟨67,-134,1⟩ | 1 | 8 | 1 | Mac | 747 | 5.3 |  |
| 268 | ⟨1,-267,1⟩ | 1 | 16 | 1 | Modal | 44851 | 9.4 |  |
| 269 | ⟨30,-90,1⟩ | 1 | 8 | ω² | Mac | 404 | 5.7 |  |
| 270 | ⟨1,-269,1⟩ | 1 | 10 | 1 | Mac | 6476 | 7.4 |  |
| 271 | ⟨17,-68,1⟩ | 1 | 4 | 1 | Mac | 375 | 5.5 |  |
| 272 | ⟨1,-271,1⟩ | 1 | 12 | ω | Mac | 4310 | 7.0 |  |
| 273 | ⟨1,-272,1⟩ | 3 | 36 | 1 | Modal | 38268 | 10.2 |  |
| 274 | ⟨11,-55,1⟩ | 1 | 3 | 1 | Mac | 405 | 11.2 |  |
| 275 | ⟨13,-35,1⟩ | 1 | 2 | ω | Mac | 368 | 5.3 |  |
| 276 | ⟨1,-275,1⟩ | 1 | 16 | 1 | Modal | 35129 | 9.0 |  |
| 277 | ⟨1,-276,1⟩ | 1 | 12 | 1 | Mac | 2713 | 7.6 |  |
| 278 | ⟨5,-19,1⟩ | 1 | 1 | ω² | Mac | 455 | 10.4 |  |
| 279 | ⟨70,-140,1⟩ | 1 | 16 | 1 | Mac | 3529 | 7.5 |  |
| 280 | ⟨1,-279,1⟩ | 1 | 10 | 1 | Modal | 15199 | 8.0 |  |
| 281 | ⟨1,-280,1⟩ | 1 | 28 | ω | Mac | 8274 | 8.2 |  |
| 282 | ⟨1,-281,1⟩ | 3 | 14 | 1 | Mac | 9607 | 7.1 |  |
| 283 | ⟨71,-142,1⟩ | 1 | 12 | 1 | Mac | 1358 | 6.8 |  |
| 284 | ⟨1,-283,1⟩ | 1 | 12 | ω | Mac | 9903 | 8.3 |  |
| 285 | ⟨1,-284,1⟩ | 1 | 16 | 1 | Mac | 6388 | 7.3 |  |
| 286 | ⟨1,-285,1⟩ | 1 | 14 | 1 | Modal | 28612 | 9.5 |  |
| 287 | ⟨2,-24,1⟩ | 1 | 3 | ω² | Mac | 353 | 5.5 |  |
| 288 | ⟨1,-17,1⟩ | 1 | 2 | 1 | Mac | 428 | 5.4 | talk or survey |
| 289 | ⟨1,-288,1⟩ | 1 | 24 | 1 | Mac | 7238 | 7.9 |  |
| 290 | ⟨1,-289,1⟩ | 1 | 12 | ω | Mac | 10722 | 7.7 |  |
| 291 | ⟨55,-74,1⟩ | 3 | 8 | 1 | Mac | 528 | 5.7 |  |
| 292 | ⟨17,-19,1⟩ | 1 | 1 | 1 | Mac | 903 | 11.7 | paper |
| 293 | ⟨6,-42,1⟩ | 1 | 4 | 1 | Mac | 454 | 6.2 |  |
| 294 | ⟨1,-293,1⟩ | 1 | 20 | 1 | Modal + Mac | 6775 | 24.2 |  |
| 295 | ⟨74,-148,1⟩ | 1 | 8 | 1 | Mac | 1226 | 6.6 |  |
| 296 | ⟨33,-99,1⟩ | 1 | 10 | ω² | Mac | 2253 | 7.6 |  |
| 297 | ⟨37,-44,1⟩ | 1 | 2 | 1 | Mac | 539 | 5.6 |  |
| 298 | ⟨1,-297,1⟩ | 1 | 12 | 1 | Mac | 7921 | 8.0 |  |
| 299 | ⟨3,-30,1⟩ | 1 | 2 | 1 | Mac | 373 | 5.8 |  |
| 300 | ⟨1,-299,1⟩ | 3 | 12 | 1 | Modal | 24952 | 9.1 |  |

## Large-d results

The table gives the companion forms of the five standard-Zauner large-d files; [README](README.md#methods) describes the constructions and representative choices.

| d | standard-Zauner form | conductor q | class number h | U_Z eigenvalue |
|---|---|---:|---:|---|
| 2208 | ⟨1,−7,1⟩ | 3 | 1 | 1 |
| 5779 | ⟨1,−3,1⟩ | 1 | 1 | 1 |
| 15128 | ⟨1,−3,1⟩ | 1 | 1 | ω² |
| 39604 | ⟨1,−3,1⟩ | 1 | 1 | 1 |
| 57966 | ⟨1,−5,1⟩ | 1 | 1 | 1 |

The separate `fibonacci/fiducial_2208.txt` has F_f symmetry, whose 16th power is type F_a; the standard 2208 file has type F_z ([Grassl–Scott, Proposition III.5](https://arxiv.org/abs/1707.02944v2)).

### Direct-ghost check

`stark_large.jl` builds the ghost directly from a double-sine series instead of an initial ghost plus Newton refinement. At d = 2208 and 5779 both routes export identical units, labels and metadata, and refine to the published files byte for byte.

### 39604 provenance

`fiducials/fiducial_39604.txt` (SHA-256 `4c20cc469b97e09b9c51995909121e3f3b1a0bf23f4d8ae61a2a1214a9e9b8b3`) was built on an Apple M2 Max (96 GiB, at most eight cores) by the prototype that this repository's large-d scripts were ported from. Its 2048-bit ghost, unit and trace checkpoints were computed before the conjugation fix; the corrected conjugation ran after it, as did reconstruction and refinement. `representatives.toml` pins its quotient generators and determinant shift. Its rebuild needs the same acceptance as the other large files: fresh construction, a byte comparison on a second platform, and both verifier modes. Until then its reproducibility from these scripts is unvalidated; the numerical checks below don't depend on it.

| stage | wall time (s) |
|---|---:|
| Ghost at 2048 bits (group, initial ghost, Newton correction; peak 79.6 GiB) | 37,578 |
| Units, eight workers | 14,734 |
| Normal bases and traces | 8,368 |
| Conjugation | 3,957 |
| Reconstruction | 9,976 |
| 640-bit refinement, eight workers | 129,408 |

`verify.py --precise` on this file took 43,942 s on one thread, with the prototype's NumPy 2.3.2 (the verifier's bytes are identical to this repository's).

### 57966 run

`fiducials/fiducial_57966.txt` (SHA-256 `262ceae8ddfdf50a5e02c3d97916f62610e5ff042598926b3a8b61ddf52eb107`) was built once, on the M2 Max (96 GiB, at most eight cores), at 2304 bits; the integer-relation step needed 1378 bits. Since then `reconstruct_large.py` has had only a behavior-preserving cleanup, checked by the byte-identical Linux rebuilds below.

| stage | dates (2026) | peak RSS (GiB) |
|---|---|---:|
| Group, ghost and units, eight threads (`--heap-size-hint=24G`) | 5–8 Nov | 14.8 |
| Normal bases, traces and conjugation | 8–9 Nov | 5.4 |
| Reconstruction | 9–10 Nov | 22.0 |
| 640-bit refinement, eight workers, 16 steps | 10–13 Nov | 1.4 (parent process) |

The laptop was in daily use, so stages were paused and restarted from checkpoints; the dates are not benchmarks. `verify.py` passed in 279 s; `verify.py --precise` took 65,068 s on one thread ([bounds](#saved-file-verification)).

### Linux rebuild measurements

`modal_rebuild.py` rebuilt the files from the current scripts in a fresh isolated copy on Linux x86, with 8 vCPU, 64 GiB reserved memory and Julia 1.13.2. Python dependencies are pinned in the script headers. Times are the runner's wall-clock seconds, excluding image building; construction includes reconstruction.

| published file | construction (s) | refinement (s) | matching SHA-256 |
|---|---:|---:|---|
| `fiducials/fiducial_2208.txt` | 1,086.42 | 269.64 | `8d922ad71a48d1c54c887d094bfa8e61f255455ddd1da4e3d4c4ff270f225d48` |
| `fibonacci/fiducial_2208.txt` | 89.81 | 285.76 | `811ca49245b4c36b0027a0ed42895073e32108104fd15880aeddf57f191eb435` |
| `fiducials/fiducial_5779.txt` | 3,760.04 | 4,979.14 | `3f5ccd8751146214272e057ed47c97201e40a3e8f481d45ad0f4ffb6cd2ac48a` |
| `fiducials/fiducial_15128.txt` | ≥ 25,895.52 | ≥ 34,932.63 | `66396d4c15e47bbac4c4d14328dcdfa7eb66586101bc839e542926a6694c66eb` |

All four rebuilt files are byte-identical to the published ones. The 15128 run took three container attempts, each resuming from checkpoints on the `sic-large-rebuild` Volume. A stage that a container restart interrupted records only its final attempt, so its times are lower bounds.

### Saved-file verification

`uv run python/verify.py` passes all 110 files. Precise bounds for the largest files:

| published file | maximum deviation, upper bound (rounded) | U_Z residual, upper bound (rounded) |
|---|---:|---:|
| `fiducials/fiducial_5779.txt` | 4.4e-153 | 9.1e-152 |
| `fiducials/fiducial_15128.txt` | 3.0e-152 | 6.6e-152 |
| `fiducials/fiducial_39604.txt` | 2.3e-152 | 5.6e-152 |
| `fiducials/fiducial_57966.txt` | 1.4e-152 | 5.4e-152 |

The 5779 and 15128 precise checks ran in one call on the M2 Max: 4,504 s wall time in total, while other jobs shared the machine. The 39604 and 57966 precise checks are in [39604 provenance](#39604-provenance) and [57966 run](#57966-run).

### Literature status

What we found for the large dimensions (sources checked 15 November 2026). Orbit equivalence with earlier SICs has not been checked.

- **2208:** reported by Grassl ([Horodecki, Rudnicki and Życzkowski, §II.1](https://arxiv.org/abs/2002.03233); [Grassl's 2021 slides](https://www.math.colostate.edu/~king/codex/slides/Grassl_2021_10_26.pdf), slide 12). We found no posted fiducial.
- **5779:** [Appleby, Bengtsson, Grassl, Harrison and McConnell](https://arxiv.org/abs/2112.05552), with [data posted by Grassl](http://s208785153.online.de/SICPOVM/5779/).
- **15128:** none found. [Bengtsson and McConnell](https://arxiv.org/html/2606.23535v1#S2) (June 2026) write that “as yet no SIC has been found.”
- **39604:** [Bengtsson, Grassl and McConnell](https://arxiv.org/abs/2403.02872). We found no posted fiducial in [Grassl's tables](http://s208785153.online.de/SICPOVM/).
- **57966:** none found.
