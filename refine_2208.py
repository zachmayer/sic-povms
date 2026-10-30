# /// script
# requires-python = "==3.13.*"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1", "python-flint==0.9.0"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Refine either d = 2208 construction to 150 significant decimal digits.

This solves the same full-space SIC equations as refine.py at 640 bits, without
imposing symmetry. At d = 2208 its dense Jacobian and ball matrices are too large:
FFT rows stream the gradient, and conjugate gradients apply the fixed float64
Gauss--Newton operator without storing it. Opposite displacement rows are paired
in the ball gradient. Only the global phase null direction is lifted. The fixed
batch order, one numerical thread, precision, tolerances and decimal formatting
preserve the two published representatives.
Unlike refine.py, retain nonzero parts below 1e-180: the original 2208 files
printed those roundoff residues, so clipping them would change their hashes.

  uv run refine_2208.py fiducial_2208_stark.txt --output fiducials/fiducial_2208.txt
  uv run refine_2208.py fiducial_2208_fibonacci_seed.txt --output fibonacci/fiducial_2208.txt
"""

import argparse
from collections.abc import Iterator
from decimal import Decimal, localcontext
import os
from pathlib import Path

# Set these before loading NumPy/SciPy, including when the caller set another count.
for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from flint import acb, arb, ctx
from scipy import fft
from scipy.sparse.linalg import LinearOperator, cg

from search import read_fiducial


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


def hp_gradient(psi: list[acb]) -> tuple[list[acb], arb, arb]:
    """J^T r and the unnormalized overlap/norm errors, with streamed FLINT DFTs."""
    d = len(psi)
    target = arb(1) / (d + 1)
    conjugates = [z.conjugate() for z in psi]
    gradient = [acb(0) for _ in psi]
    worst = arb(0)
    norm = sum((z.real * z.real + z.imag * z.imag for z in psi), arb(0))
    for a in range(d // 2 + 1):
        weight = 1 if a in (0, d // 2) else 2
        indices = [(j + a) % d for j in range(d)]
        c = acb.dft([conjugates[k] * z for k, z in zip(indices, psi, strict=True)])
        residual = [z.real * z.real + z.imag * z.imag - target for z in c]
        if a == 0:
            residual[0] = c[0].real * c[0].real - 1
        worst = max(worst, max(abs(z).mid() for z in residual[1 if a == 0 else 0:]))
        h = acb.dft([2 * d * weight * r * z for r, z in zip(residual, c, strict=True)], inverse=True)
        for j, k in enumerate(indices):
            gradient[j] += h[j] * psi[k]
            gradient[k] += h[j].conjugate() * psi[j]
    return gradient, worst, abs(norm - 1).mid()


def scott(x: arb, digits: int) -> str:
    """Round a ball midpoint to the original Scott--Grassl decimal convention."""
    mantissa, exponent = x.mid().man_exp()
    with localcontext() as context:
        context.prec = digits + 30
        value = Decimal(int(mantissa)) * Decimal(2) ** int(exponent)
        if value == 0:
            return "+." + "0" * digits + "e0"
        a, b = format(abs(value), f".{digits - 1}e").split("e")
        return ("+" if value > 0 else "-") + "." + a.replace(".", "") + "e" + str(int(b) + 1)


def write_hp(path: Path, psi: list[acb]) -> None:
    """Normalize and fix the first component positive before rounding to 150 digits."""
    norm = sum((z.real * z.real + z.imag * z.imag for z in psi), arb(0)).sqrt()
    phase = psi[0].conjugate() / abs(psi[0])
    normalized = [z * phase / norm for z in psi]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(scott(z.real, 150) + scott(arb(0) if j == 0 else z.imag, 150) + "i\n"
                            for j, z in enumerate(normalized)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("fiducial", type=Path, help="the 17-digit construction seed")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ctx.prec = 640
    ctx.threads = 1
    strings = read_fiducial(args.fiducial)
    if len(strings) != 2208:
        raise ValueError("This refinement driver is for d = 2208")
    psi = [acb(arb(re), arb(im)).mid() for re, im in strings]
    p0 = np.array([complex(float(re), float(im)) for re, im in strings])
    d = len(psi)
    phase = np.concatenate([-p0.imag, p0.real])
    phase /= np.linalg.norm(phase)

    def normal(x: np.ndarray) -> np.ndarray:
        z = normal_product(p0, x[:d] + 1j * x[d:], batch=32)
        return np.concatenate([z.real, z.imag]) + phase * np.dot(phase, x)

    operator = LinearOperator((2 * d, 2 * d), matvec=normal, dtype=np.float64)
    for step in range(22):
        gradient, error, norm_error = hp_gradient(psi)
        print(f"step {step}: overlap error {float(error):.1e}, norm error {float(norm_error):.1e}", flush=True)
        if max(float(error), float(norm_error)) < 1e-180:
            write_hp(args.output, psi)
            print(f"wrote {args.output}")
            return
        rhs = -np.array([float(z.real.mid()) for z in gradient] + [float(z.imag.mid()) for z in gradient])
        scale = np.max(abs(rhs))
        rhs /= scale
        delta, info = cg(operator, rhs, rtol=1e-12, atol=0, maxiter=600)
        if info:
            raise RuntimeError(f"Normal solve failed: {info}")
        delta *= scale
        psi = [(z + acb(delta[j], delta[d + j])).mid() for j, z in enumerate(psi)]
    raise RuntimeError("Refinement did not converge in 22 steps")


if __name__ == "__main__":
    main()
