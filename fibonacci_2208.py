# /// script
# requires-python = "==3.13.*"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Rebuild the d=2208 Fibonacci-family double-precision seed, seed 2348, trial 3.

This F_a representative needs an alignment-constrained search rather than stark.jl's
standard-Zauner reconstruction. CRT factors of the Fibonacci antiunitary give 96 real
coordinates; both even-dimensional alignment relations give 118 quadratic constraints.
Projected L-BFGS with least-squares retraction minimizes the full SIC loss. FFT rows
replace the dense overlap arrays of search.py/refine.py at this dimension. No cached
basis, checkpoint, or previous trial is read. All arithmetic choices and the local RNG
are fixed; BLAS and FFT calls use one thread. refine_2208.py supplies the 150 digits.

fiducial_48g.txt is Grassl--Scott's published 48g, unchanged numerical text from
https://arxiv.org/src/1707.02944 (SHA-256 below). The alignment relations are Eqs.15--16
of https://arxiv.org/abs/1905.09737. The antiunitary is Eqs.24--25 of 1707.02944.

  uv run fibonacci_2208.py
"""

import argparse
import hashlib
import os
from pathlib import Path

# These must precede numerical imports: the seeded optimizer's representative depends
# on LAPACK's rounding, so thread-count changes must not change its arithmetic path.
for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy import fft
from scipy.linalg import eigh, schur, svd
from scipy.optimize import least_squares

from search import read_fiducial

D = 2208
K = 8
INPUT_SHA256 = "f0b3517ccf6e9ebe40eea6cee62b436f80226e5bf99280a2f0b71f894dc9d090"


def antiunitary(psi: np.ndarray) -> np.ndarray:
    """A psi = tau^(r^2) sqrt(d) ifft(conj(psi)), including even-d label signs."""
    r = np.arange(D, dtype=np.int64)
    chirp = np.exp(1j * np.pi * ((D + 1) * r * r % (2 * D)) / D)
    return chirp * fft.ifft(np.conj(psi), axis=0, workers=1) * np.sqrt(D)


def tensor_apply(matrices: list[np.ndarray], vector: np.ndarray) -> np.ndarray:
    for axis, matrix in enumerate(matrices):
        vector = np.moveaxis(np.tensordot(matrix, vector, axes=(1, axis)), 0, axis)
    return vector


def fixed_basis() -> np.ndarray:
    """Real orthonormal basis for A psi=psi, using local Schur factors of A^2."""
    factors = [32, 3, 23]
    projectors, local_bases = [], []
    eigenvalues = np.ones(factors, dtype=complex)
    for axis, n in enumerate(factors):
        m = D // n
        t = pow(m, -1, n)
        e = m * t
        r = np.arange(n, dtype=np.int64)
        j = e * r % D
        chirp = np.exp(1j * np.pi * ((D + 1) * j * j % (2 * D)) / D)
        fourier = np.exp(2j * np.pi * ((t * r[:, None] * r) % n) / n) / np.sqrt(n)
        matrix = chirp[:, None] * fourier
        triangular, q = schur(matrix @ matrix.conj(), output="complex")
        if np.max(np.abs(triangular - np.diag(np.diag(triangular)))) > 1e-10:
            raise ValueError("Local unitary did not diagonalize")
        shape = [1] * len(factors)
        shape[axis] = n
        eigenvalues *= np.diag(triangular).reshape(shape)
        local_bases.append(q)
        projectors.append(e)
    grid = np.indices(factors, dtype=np.int64)
    indices = sum(e * grid[i] for i, e in enumerate(projectors)) % D
    keep = np.abs(eigenvalues - 1) < 1e-8
    rank = int(np.count_nonzero(keep))
    if rank != 96:
        raise ValueError(f"Expected fixed-space rank 96, got {rank}")

    def decode(c: np.ndarray) -> np.ndarray:
        coefficients = np.zeros(factors, dtype=complex)
        coefficients[keep] = c
        vector = np.empty(D, dtype=complex)
        vector[indices] = tensor_apply(local_bases, coefficients)
        return (vector + antiunitary(vector)) / 2

    columns = [decode(z) for z in np.eye(rank)]
    columns += [decode(1j * z) for z in np.eye(rank)]
    b = np.array(columns).T
    u, singular, _ = svd(np.concatenate([b.real, b.imag]), full_matrices=False, lapack_driver="gesdd")
    if np.sum(singular > 1e-8) != rank:
        raise ValueError("Real fixed-space rank disagrees with local Schur factors")
    return u[:D, :rank] + 1j * u[D:, :rank]


def alignment_targets(source: Path) -> dict[tuple[int, int], complex]:
    """Both even-d alignment relations, with the identity map on 48g's labels."""
    if hashlib.sha256(source.read_bytes()).hexdigest() != INPUT_SHA256:
        raise ValueError(f"Unexpected Grassl--Scott 48g input: {source}")
    psi = np.array([complex(float(re), float(im)) for re, im in read_fiducial(source)])
    s = 48
    a, b = np.indices((s, s))
    c = fft.ifft(psi * psi[(a + b) % s].conj(), axis=1, workers=1) * s
    c *= np.exp(1j * np.pi * ((s + 1) * a * b % (2 * s)) / s)
    small = c * np.sqrt(s + 1)
    small[0, 0] = 1
    values = (-1.0) ** ((a + 1) * (b + 1)) * small[a, b] ** 2 / np.sqrt(D + 1)
    values[0, 0] = 1
    targets = {(int((s - 2) * p), int((s - 2) * q)): values[p, q]
               for p, q in zip(a.ravel(), b.ravel(), strict=True)}
    for p in range(s - 2):
        for q in range(s - 2):
            value = 1 if (p, q) == (0, 0) else -(-1.0) ** ((p + 1) * (q + 1)) / np.sqrt(D + 1)
            key = (s * p, s * q)
            if key in targets and abs(targets[key] - value) >= 1e-12:
                raise ValueError("The alignment relations disagree")
            targets[key] = value
    return targets


class Alignment:
    """Independent orthonormal quadratic polynomials for the alignment relations."""

    def __init__(self, source: Path) -> None:
        targets = alignment_targets(source)
        self.b = fixed_basis()
        unseen, representatives = set(targets), []
        while unseen:
            point = min(unseen)
            representatives.append(point)
            a, b = point
            for _ in range(6 * K):
                unseen.discard((a, b))
                unseen.discard((-a % D, -b % D))
                a, b = b, (a + b) % D
            if (a, b) != point:
                raise ValueError("Fibonacci orbit did not close")
        operators, values = [], []
        j = np.arange(D)
        for a, b in representatives:
            phase = np.exp(1j * np.pi * ((2 * b * j + (D + 1) * a * b) % (2 * D)) / D)
            op = self.b[(j + a) % D].conj().T @ (phase[:, None] * self.b)
            op = (op + op.T) / 2
            scale = 1 if (a, b) == (0, 0) else np.sqrt(D + 1)
            operators.extend([scale * op.real, scale * op.imag])
            values.extend([scale * targets[a, b].real, scale * targets[a, b].imag])
        operators, values = np.array(operators), np.array(values)
        flat = operators.reshape(len(operators), -1)
        eigenvalues, u = eigh(flat @ flat.T)
        keep = eigenvalues > eigenvalues[-1] * 1e-12
        if np.max(abs(u[:, ~keep].T @ values), initial=0) > 1e-7:
            raise ValueError("Inconsistent alignment polynomial relations")
        transform = (u[:, keep] / np.sqrt(eigenvalues[keep])).T
        self.q = (transform @ flat).reshape(-1, self.b.shape[1], self.b.shape[1])
        self.t = transform @ values

    def residual(self, x: np.ndarray) -> np.ndarray:
        return (self.q @ x) @ x - self.t

    def jacobian(self, x: np.ndarray) -> np.ndarray:
        return 2 * (self.q @ x)

    def tangent(self, x: np.ndarray) -> np.ndarray:
        _, singular, v = svd(self.jacobian(x), full_matrices=False, lapack_driver="gesdd")
        return v[singular < 1e-7].T

    def retract(self, x: np.ndarray) -> np.ndarray:
        fit = least_squares(self.residual, x, jac=self.jacobian, method="lm", max_nfev=100,
                            ftol=1e-13, gtol=1e-13, xtol=1e-13)
        if np.max(abs(self.residual(fit.x))) > 1e-9:
            raise RuntimeError("Alignment retraction failed")
        return fit.x


def objective(b: np.ndarray, x: np.ndarray) -> tuple[float, np.ndarray]:
    """Half the FFT rows give full SIC loss, since |C_(a,b)|=|C_(-a,-b)|."""
    psi = b @ x
    shifts = np.arange(D // 2 + 1)
    weights = np.full(len(shifts), 2.0)
    weights[0] = weights[-1] = 1
    loss, gradient = 0.0, np.zeros_like(psi)
    k = np.arange(D, dtype=np.int64)
    for start in range(0, shifts.size, 32):
        a = shifts[start:start + 32]
        w = weights[start:start + 32, None]
        indices = (a[:, None] + k) % D
        c = fft.fft(np.conj(psi[indices]) * psi, axis=-1, workers=1)
        residual = np.abs(c) ** 2 - 1 / (D + 1)
        if a[0] == 0:
            residual[0, 0] = abs(c[0, 0]) ** 2 - 1
        loss += float(np.sum(w * residual * residual)) / 2
        h = D * fft.ifft(2 * (w * residual) * c, axis=-1, workers=1)
        g = np.sum(h * psi[indices], axis=0)
        for row, shift in enumerate(a):
            g += np.roll(np.conj(h[row]) * psi, int(shift))
        gradient += g
    return loss, (b.conj().T @ gradient).real


def minimize(model: Alignment, start: np.ndarray) -> np.ndarray:
    """Projected L-BFGS, retaining the original trial's operation order."""
    x = model.retract(start)
    loss, gradient = objective(model.b, x)
    tangent = model.tangent(x)
    projected = tangent @ (tangent.T @ gradient)
    history = []
    for iteration in range(500):
        q, alphas, pairs = projected.copy(), [], []
        for s, y in history:
            s, y = tangent @ (tangent.T @ s), tangent @ (tangent.T @ y)
            sy = np.dot(s, y)
            if sy > 1e-14 * np.linalg.norm(s) * np.linalg.norm(y):
                pairs.append((s, y, 1 / sy))
        for s, y, rho in reversed(pairs):
            alpha = rho * np.dot(s, q)
            alphas.append(alpha)
            q -= alpha * y
        gamma = np.dot(pairs[-1][0], pairs[-1][1]) / np.dot(pairs[-1][1], pairs[-1][1]) if pairs else 1
        direction = gamma * q
        for (s, y, rho), alpha in zip(pairs, reversed(alphas), strict=True):
            direction += s * (alpha - rho * np.dot(y, direction))
        direction = -direction
        slope = np.dot(projected, direction)
        if slope >= 0:
            direction, slope, history = -projected, -np.dot(projected, projected), []
        if np.linalg.norm(projected) < 1e-11 or loss < 1e-22:
            print(f"trial 3: {iteration} iterations, loss {loss:.3e}", flush=True)
            return x
        step = min(1.0, 0.5 / max(np.linalg.norm(direction), 1e-15))
        for _ in range(24):
            newx = model.retract(x + step * direction)
            newloss, newgradient = objective(model.b, newx)
            if newloss <= loss + 1e-4 * step * slope:
                break
            step *= 0.5
        else:
            raise RuntimeError("Projected L-BFGS line search failed")
        newtangent = model.tangent(newx)
        newprojected = newtangent @ (newtangent.T @ newgradient)
        history.append((newx - x, newprojected - newtangent @ (newtangent.T @ projected)))
        history = history[-15:]
        x, loss, tangent, projected = newx, newloss, newtangent, newprojected
    raise RuntimeError("No convergence in 500 iterations")


def write_seed(path: Path, psi: np.ndarray) -> None:
    psi = psi / np.linalg.norm(psi)
    psi *= np.conj(psi[0]) / abs(psi[0])
    psi[0] = psi[0].real

    def number(x: float) -> str:
        if x == 0:
            return "+.00000000000000000e0"
        mantissa, exponent = f"{abs(x):.16e}".split("e")
        return f"{'+' if x > 0 else '-'}.{mantissa.replace('.', '')}e{int(exponent) + 1}"

    path.write_text("".join(number(z.real) + number(z.imag) + "i\n" for z in psi))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, default=Path("fiducial_2208_fibonacci_seed.txt"))
    args = parser.parse_args()
    model = Alignment(Path(__file__).with_name("fiducial_48g.txt"))
    rng = np.random.default_rng(2348)
    for _ in range(4):
        start = rng.normal(size=model.b.shape[1])
        start /= np.linalg.norm(start)
    x = minimize(model, start)
    psi = model.b @ x
    # A failed search must never silently become refinement input.
    loss, _ = objective(model.b, x)
    if loss >= 1e-20:
        raise RuntimeError(f"Trial 3 did not reach a near-SIC: loss {loss}")
    write_seed(args.output, psi)
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
