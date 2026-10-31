# /// script
# requires-python = "==3.13.*"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1", "python-flint==0.9.0"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""640-bit numerical SIC refinement with automatic atomic restart checkpoints.

FLINT FFTs stream over one representative of the conjugate displacement pairs.
The double-precision normal operator stays anchored to the ORIGINAL input file;
only its phase null direction is lifted. The input must already be near a SIC.
No double-precision invariant basis constrains the high-precision answer.

Checkpoints store exact binary midpoints, never normalized decimal exports.
They retain gradient row progress and CG approximations. With --workers 2..8,
independent gradient shards run in separate processes and checkpoint separately.
The normal solve is serial. All BLAS/FFT thread counts must be one per process
(this script sets these). Final files require an independent all-overlap and
symmetry check (standard Zauner or Fibonacci, according to the output folder).
This is numerical refinement. All published large-d inputs select representatives
from an unpublished prototype run via representatives.toml, retained so its Mac
runs and Linux rebuilds compare byte for byte.
"""

import argparse
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import re
import resource
import sys
import time

for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from flint import acb, arb, ctx
from scipy import fft
from scipy.sparse.linalg import LinearOperator, cg

from search import fiducial_text, read_fiducial


def batches(d: int, batch: int) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Displacement shifts and their component indices, in a fixed order."""
    k = np.arange(d, dtype=np.int64)
    for start in range(0, d, batch):
        a = np.arange(start, min(start + batch, d), dtype=np.int64)
        yield a, (a[:, None] + k) % d


def normal_product(psi: np.ndarray, tangent: np.ndarray, batch: int = 32) -> np.ndarray:
    """J^T J tangent for |C_ab|^2, packed as Re g + i Im g, using O(batch*d) memory."""
    gradient = np.zeros_like(psi)
    for a, indices in batches(psi.size, batch):
        c = fft.fft(np.conj(psi[indices]) * psi, axis=-1, workers=1)
        dc = fft.fft(np.conj(tangent[indices]) * psi + np.conj(psi[indices]) * tangent,
                     axis=-1, workers=1)
        dr = 2 * np.real(np.conj(c) * dc)
        h = psi.size * fft.ifft(2 * dr * c, axis=-1, workers=1)
        batch_gradient = np.sum(h * psi[indices], axis=0)
        for row, shift in enumerate(a):
            batch_gradient += np.roll(np.conj(h[row]) * psi, int(shift))
        gradient += batch_gradient
    return gradient



ROW = re.compile(r"([+-]\.\d+e-?\d+)([+-]\.\d+e-?\d+)i")


def read_hp(path):
    answer = []
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip():
            continue
        match = ROW.fullmatch(line.strip())
        if match is None:
            raise ValueError(f"Invalid Scott--Grassl row {number}: {path}")
        a, b = match.groups()
        answer.append(acb(arb(a.replace(".", "0.", 1)),
                          arb(b.replace(".", "0.", 1))).mid())
    if len(answer) < 2:
        raise ValueError("At least two components required")
    return answer


def atomic_text(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".pending")
    temporary.write_text(text)
    temporary.replace(path)


def fraction(x: arb) -> Fraction:
    """Use the same zero threshold as refine.py, then retain the binary midpoint."""
    mantissa, exponent = x.mid().man_exp()
    value = Fraction(int(mantissa)) * Fraction(2) ** int(exponent)
    return Fraction(0) if abs(value) < Fraction(1, 10**180) else value


def write_hp(path: Path, psi: list[acb]) -> None:
    """Normalize and fix the first component positive before rounding to 150 digits."""
    norm = sum((z.real * z.real + z.imag * z.imag for z in psi), arb(0)).sqrt()
    phase = psi[0].conjugate() / abs(psi[0])
    normalized = [z * phase / norm for z in psi]
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [(fraction(z.real), Fraction(0) if j == 0 else fraction(z.imag))
            for j, z in enumerate(normalized)]
    temporary = path.with_suffix(path.suffix + ".pending")
    temporary.write_bytes(fiducial_text(rows, 150).encode("ascii"))
    temporary.replace(path)



def pack_real(x):
    mantissa, exponent = x.mid().man_exp()
    return [str(mantissa), int(exponent)]


def unpack_real(pair):
    return arb((int(pair[0]), pair[1]))


def pack_vector(vector):
    return [[pack_real(z.real), pack_real(z.imag)] for z in vector]


def unpack_vector(vector):
    return [acb(unpack_real(a), unpack_real(b)) for a, b in vector]


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def save_state(path, state):
    atomic_text(path, json.dumps(state, separators=(",", ":"), allow_nan=False) + "\n")


def peak_rss_bytes():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(rss if sys.platform == "darwin" else rss*1024)


def hp_gradient(psi, rows=None, *, start_row=0, gradient=None, worst=None,
                checkpoint=None, checkpoint_rows=64):
    """J^T r over [start_row, rows), or complete with a saved prefix.

    Row a=d//2 has weight one only for EVEN d. For odd d it is paired with
    d-a and therefore has weight two. Midpoints are taken after each row,
    identically with/without checkpoints. The independent checker uses bounds.
    """
    d = len(psi)
    stop = d//2 + 1 if rows is None else min(rows, d//2 + 1)
    if not 0 <= start_row <= stop:
        raise ValueError("Invalid gradient row interval")
    target = arb(1)/(d+1)
    conjugates = [z.conjugate() for z in psi]
    gradient = [acb(0) for _ in psi] if gradient is None else gradient
    worst = arb(0) if worst is None else worst
    norm = sum((z.real*z.real + z.imag*z.imag for z in psi), arb(0))
    for a in range(start_row, stop):
        weight = 1 if a == 0 or (d % 2 == 0 and a == d//2) else 2
        indices = [(j+a) % d for j in range(d)]
        overlaps = acb.dft([conjugates[k]*z for k, z in zip(indices, psi, strict=True)])
        residual = [z.real*z.real + z.imag*z.imag - target for z in overlaps]
        if a == 0:
            residual[0] = overlaps[0].real*overlaps[0].real - 1
        worst = max(worst, max(abs(z).mid() for z in residual[1 if a == 0 else 0:]))
        h = acb.dft([2*d*weight*r*z for r, z in zip(residual, overlaps, strict=True)],
                    inverse=True)
        for j, k in enumerate(indices):
            gradient[j] += h[j]*psi[k]
            gradient[k] += h[j].conjugate()*psi[j]
        gradient = [z.mid() for z in gradient]
        if checkpoint is not None and ((a+1) % checkpoint_rows == 0 or a+1 == stop):
            checkpoint(a+1, gradient, worst)
    return gradient, worst, abs(norm-1).mid()


def gradient_shard(packed_psi, identity, first, stop, path, checkpoint_rows):
    """Independently callable, resumable worker; no network or global state."""
    ctx.prec = identity["bits"]
    ctx.threads = 1
    path = Path(path)
    shard_identity = {**identity, "psi_sha256": digest(packed_psi), "first": first, "stop": stop}
    if path.exists():
        state = json.loads(path.read_text())
        if state["identity"] != shard_identity:
            raise ValueError(f"Gradient shard input/settings mismatch: {path}")
        gradient = unpack_vector(state["gradient"])
    else:
        gradient = [acb(0) for _ in packed_psi]
        state = {"identity": shard_identity, "next_row": first,
                 "gradient": pack_vector(gradient), "worst": pack_real(arb(0)), "seconds": 0.0}
        save_state(path, state)
    if state["next_row"] == stop:
        return state
    start = time.perf_counter()
    prior_seconds = state["seconds"]

    def checkpoint(next_row, partial, worst):
        state.update(next_row=next_row, gradient=pack_vector(partial), worst=pack_real(worst),
                     seconds=prior_seconds + time.perf_counter()-start)
        save_state(path, state)

    hp_gradient(unpack_vector(packed_psi), stop, start_row=state["next_row"], gradient=gradient,
                worst=unpack_real(state["worst"]), checkpoint=checkpoint,
                checkpoint_rows=checkpoint_rows)
    return state


def merge_gradient_shards(parts, identity, packed_psi):
    """Reject gaps, duplicate rows, stale iterates and incomplete shards."""
    parts = sorted(parts, key=lambda part: part["identity"]["first"])
    cursor = 0
    gradient = [acb(0) for _ in packed_psi]
    worst = arb(0)
    psi_hash = digest(packed_psi)
    for part in parts:
        key = part["identity"]
        if (any(key.get(k) != v for k, v in identity.items())
                or key["psi_sha256"] != psi_hash or key["first"] != cursor
                or key["stop"] <= cursor or part["next_row"] != key["stop"]):
            raise ValueError("Invalid, overlapping, incomplete or incompatible gradient shards")
        gradient = [(a+b).mid() for a, b in zip(gradient, unpack_vector(part["gradient"]), strict=True)]
        worst = max(worst, unpack_real(part["worst"]))
        cursor = key["stop"]
    if cursor != len(packed_psi)//2 + 1:
        raise ValueError("Gradient shards do not cover every representative row")
    return gradient, worst


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=22, help="Total iteration limit, including resumed ones")
    parser.add_argument("--bits", type=int, default=640)
    parser.add_argument("--workers", type=int, default=1, help="Gradient workers (1..8); CG stays single-core")
    parser.add_argument("--checkpoint", type=Path)
    parser.add_argument("--checkpoint-rows", type=int, default=64)
    parser.add_argument("--checkpoint-cg-iterations", type=int, default=16)
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--benchmark-rows", type=int, default=16)
    args = parser.parse_args()
    if args.bits < 640 or not 1 <= args.workers <= 8 or min(
            args.steps, args.checkpoint_rows, args.checkpoint_cg_iterations, args.benchmark_rows) < 1:
        parser.error("Use at least 640 bits, 1..8 workers and positive iteration/row counts")
    if args.file.resolve() == args.output.resolve():
        parser.error("Input and output must differ; input anchors the normal operator")
    ctx.prec = args.bits
    ctx.threads = 1
    psi = read_hp(args.file)
    d = len(psi)
    if args.benchmark:
        start = time.perf_counter()
        _, worst, norm = hp_gradient(psi, args.benchmark_rows)
        print(json.dumps({"rows": min(args.benchmark_rows, d//2+1), "d": d,
                          "bits": args.bits, "seconds": time.perf_counter()-start,
                          "error": float(worst), "norm_error": float(norm),
                          "peak_rss_bytes": peak_rss_bytes()}), flush=True)
        return
    checkpoint_path = args.checkpoint or args.output.with_suffix(".checkpoint.json")
    if checkpoint_path.resolve() in (args.file.resolve(), args.output.resolve()):
        parser.error("Checkpoint must differ from input and output")
    identity = {"input_sha256": hashlib.sha256(args.file.read_bytes()).hexdigest(),
                "bits": args.bits, "d": d, "workers": min(args.workers, d//2+1),
                "code_sha256": hashlib.sha256(Path(__file__).read_bytes()
                                               + Path(__file__).with_name("search.py").read_bytes()).hexdigest()}
    if checkpoint_path.exists():
        state = json.loads(checkpoint_path.read_text())
        if state["identity"] != identity:
            raise ValueError("Checkpoint input, precision, workers or implementation changed; use another checkpoint")
        psi = unpack_vector(state["psi"])
        print(json.dumps({"event": "resume", "step": state["step"], "phase": state["phase"],
                          "next_row": state["next_row"], "checkpoint": str(checkpoint_path)}), flush=True)
    else:
        state = {"identity": identity, "step": 0, "phase": "gradient",
                 "psi": pack_vector(psi), "gradient": None, "next_row": 0,
                 "worst": pack_real(arb(0)), "gradient_seconds": 0.0, "cg_x": None,
                 "normal_seconds": 0.0, "normal_calls": 0}
        save_state(checkpoint_path, state)
    if state["phase"] == "done":
        write_hp(args.output, psi)
        print(json.dumps({"event": "already_converged", "output": str(args.output)}), flush=True)
        return
    p0 = np.array([complex(float(re), float(im)) for re, im in read_fiducial(args.file)])
    phase = np.concatenate([-p0.imag, p0.real])
    phase /= np.linalg.norm(phase)
    calls = 0

    def normal(x):
        nonlocal calls
        calls += 1
        z = normal_product(p0, x[:d] + 1j*x[d:], batch=32)
        return np.concatenate([z.real, z.imag]) + phase*np.dot(phase, x)

    operator = LinearOperator((2*d, 2*d), matvec=normal, dtype=np.float64)
    log = args.output.with_suffix(".jsonl")
    for step in range(state["step"], args.steps):
        gradient = None if state["gradient"] is None else unpack_vector(state["gradient"])
        if state["phase"] == "gradient":
            start = time.perf_counter()
            prior_seconds = state["gradient_seconds"]
            if identity["workers"] == 1:

                def checkpoint_gradient(next_row, partial, worst):
                    state.update(next_row=next_row, gradient=pack_vector(partial), worst=pack_real(worst),
                                 gradient_seconds=prior_seconds + time.perf_counter()-start)
                    save_state(checkpoint_path, state)

                gradient, error, norm_error = hp_gradient(
                    psi, start_row=state["next_row"], gradient=gradient,
                    worst=unpack_real(state["worst"]), checkpoint=checkpoint_gradient,
                    checkpoint_rows=args.checkpoint_rows)
            else:
                count = identity["workers"]
                rows = d//2 + 1
                packed = state["psi"]
                context = multiprocessing.get_context("spawn")
                with ProcessPoolExecutor(max_workers=count, mp_context=context) as pool:
                    jobs = [pool.submit(gradient_shard, packed, identity, rows*k//count, rows*(k+1)//count,
                                        str(checkpoint_path.with_name(checkpoint_path.name
                                            + f".step{step}.part{k}.json")), args.checkpoint_rows)
                            for k in range(count)]
                    parts = [job.result() for job in jobs]
                gradient, error = merge_gradient_shards(parts, identity, packed)
                norm_error = abs(sum((z.real*z.real + z.imag*z.imag for z in psi), arb(0))-1).mid()
                state["gradient_core_seconds"] = sum(part["seconds"] for part in parts)
            state.update(phase="solve", gradient=pack_vector(gradient), worst=pack_real(error),
                         norm_error=pack_real(norm_error), next_row=d//2+1,
                         gradient_seconds=prior_seconds + time.perf_counter()-start)
            save_state(checkpoint_path, state)
        error = unpack_real(state["worst"])
        norm_error = unpack_real(state["norm_error"])
        row = {"step": step, "bits": ctx.prec, "max_error": float(error),
               "norm_error": float(norm_error), "gradient_seconds": state["gradient_seconds"],
               "peak_rss_bytes": peak_rss_bytes(), "workers": identity["workers"]}
        print(json.dumps(row), flush=True)
        with log.open("a") as stream:
            stream.write(json.dumps(row) + "\n")
        if max(float(error), float(norm_error)) < 1e-180:
            write_hp(args.output, psi)
            state["phase"] = "done"
            save_state(checkpoint_path, state)
            print(json.dumps({"event": "converged", "output": str(args.output),
                              "step": step, "digits": 150}), flush=True)
            return
        rhs = -np.array([float(z.real.mid()) for z in gradient]
                        + [float(z.imag.mid()) for z in gradient])
        scale = np.max(abs(rhs))
        if not np.isfinite(scale) or scale == 0:
            raise RuntimeError("Gradient is nonfinite or vanished before convergence")
        rhs /= scale
        calls = 0
        cg_iterations = 0
        start = time.perf_counter()
        prior_seconds = state["normal_seconds"]
        prior_calls = state["normal_calls"]

        def checkpoint_cg(x):
            nonlocal cg_iterations
            cg_iterations += 1
            if cg_iterations % args.checkpoint_cg_iterations == 0:
                state.update(cg_x=x.tolist(), normal_seconds=prior_seconds + time.perf_counter()-start,
                             normal_calls=prior_calls + calls)
                save_state(checkpoint_path, state)

        x0 = None if state["cg_x"] is None else np.array(state["cg_x"])
        delta, info = cg(operator, rhs, x0=x0, rtol=1e-12, atol=0, maxiter=600,
                         callback=checkpoint_cg)
        if info:
            raise RuntimeError(f"Normal solve failed: {info}; last CG checkpoint retained")
        delta *= scale
        psi = [(z + acb(delta[j], delta[d+j])).mid() for j, z in enumerate(psi)]
        print(json.dumps({"step": step, "normal_calls": prior_calls + calls,
                          "normal_seconds": prior_seconds + time.perf_counter()-start}), flush=True)
        state.update(step=step+1, phase="gradient", psi=pack_vector(psi), gradient=None,
                     next_row=0, worst=pack_real(arb(0)), gradient_seconds=0.0,
                     cg_x=None, normal_seconds=0.0, normal_calls=0)
        save_state(checkpoint_path, state)
    raise RuntimeError("Refinement did not converge within --steps; checkpoint retained")


if __name__ == "__main__":
    main()
