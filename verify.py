# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.4", "python-flint==0.9.0"]
# [tool.uv]
# exclude-newer = "2026-10-30T00:00:00Z"
# ///
"""Check every fiducial_<d>.txt in fiducials/ and fibonacci/, or the given files.

A unit vector psi in C^d is a SIC fiducial iff |<psi|X^a Z^b|psi>|^2 = 1/(d+1)
for every (a,b) != (0,0), where X|k> = |k+1 mod d> and
Z|k> = exp(2 pi i k/d)|k>. We also check the filename's d and psi_0 real positive.
Scott--Grassl's conventions (arXiv:0910.5784, arXiv:1703.03993) require
U_Z[r,s] = exp(i pi (d-1)/12)/sqrt(d) exp(i pi (2rs + (d+1)r^2)/d)
to have psi as an eigenvector (checked against their d = 52,55,60,67,70,88).

For files in fibonacci/, check Grassl--Scott's antiunitary F_f = [[0,1],[1,1]]
(arXiv:1707.02944v2, Eqs. 24--25) instead:
A psi[r] = sum_s tau^(r^2 + 2rs) conjugate(psi[s])/sqrt(d), tau = -exp(pi i/d).
Its eigenphase is arbitrary. At d = L_(2k)+1, F_f^(2k) has order 3 and is type
F_a when 4 divides k (Proposition III.5); d = 2208 has k = 8.
Both chirps use (d+1)r^2 modulo 2d, retaining the signs for even d.

The negative DFT of conjugate(psi[k+a])*psi[k] gives the overlaps at b -> -b.
One row at a time checks all d^2 overlaps in O(d^2 log d) time, O(d) space.
The symmetry uses a positive DFT, obtained by reversing the negative DFT.
Precise transforms use Bluestein convolution with cached chirps and power-of-two
FLINT FFTs, so composite and prime dimensions follow the same algorithm.
This independent judge imports no other repository scripts.

uv run verify.py [FILE ...]            # NumPy double precision; threshold 1e-12
uv run verify.py [FILE ...] --precise  # FLINT 640-bit balls; threshold 1e-147
In --precise, max deviation (including norm) and symmetry residual are rigorous
upper bounds for the decimal input; PASS requires both bounds below threshold.
The printed bounds are rounded to two significant digits, as in default mode.
"""

import cmath
import re
import sys
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx


def squared_abs(z):
    """Squared modulus without a square root or loss of ball enclosure."""
    return z.real * z.real + z.imag * z.imag


def fourier(d, precise):
    """Unnormalized negative DFT; cache one Bluestein kernel per precise file."""
    if not precise:
        return np.fft.fft
    n = 1 << (2 * d - 2).bit_length()  # At least 2d-1; only the middle d convolution terms are used.
    chirp = [(acb(-(k * k % (2 * d))) / d).exp_pi_i() for k in range(d)]
    kernel = [chirp[abs(k)].conjugate() for k in range(1 - d, d)]
    kernel_fft = acb.dft(kernel + [acb(0)] * (n - len(kernel)))

    def transform(vector):
        # exp(-2 pi i jk/d) = chirp[j] chirp[k] conjugate(chirp[j-k]).
        forward = acb.dft([v * c for v, c in zip(vector, chirp, strict=True)] + [acb(0)] * (n - d))
        convolution = acb.dft([v * c for v, c in zip(forward, kernel_fft, strict=True)], inverse=True)
        return [convolution[k + d - 1] * chirp[k] for k in range(d)]

    return transform


def check(psi, precise, fft):
    """Return the maximum SIC overlap and unit-norm error."""
    d = len(psi)
    if precise:
        target = arb(1) / (d + 1)
        worst = abs(sum(squared_abs(z) for z in psi) - 1).upper()
        conjugate = [z.conjugate() for z in psi]
    else:
        target = 1.0 / (d + 1)
        worst = abs(sum(abs(z) ** 2 for z in psi.tolist()) - 1)
        conjugate = psi.conjugate()
    for a in range(d):
        if precise:
            spectrum = fft([conjugate[(k + a) % d] * psi[k] for k in range(d)])
            errors = (abs(squared_abs(z) - target).upper() for z in spectrum[1 if a == 0 else 0:])
            worst = max(worst, max(errors, default=worst))
        else:
            spectrum = fft(np.roll(conjugate, -a) * psi)
            errors = np.abs(np.abs(spectrum[1 if a == 0 else 0:]) ** 2 - target)
            worst = max(worst, float(np.max(errors, initial=0)))
    return worst


def symmetry(psi, fibonacci, precise, fft):
    """Apply the stated unitary/antiunitary; all phase labels live modulo 2d."""
    d = len(psi)
    exponents = [((d + 1) * r * r) % (2 * d) for r in range(d)]
    if precise:
        vector = [z.conjugate() for z in psi] if fibonacci else psi
        transformed = fft(vector)
        scale = (acb(1) if fibonacci else (acb(d - 1) / 12).exp_pi_i()) / arb(d).sqrt()
        image = [scale * (acb(e) / d).exp_pi_i() * transformed[-r % d]
                 for r, e in enumerate(exponents)]
        lam = sum(p.conjugate() * q for p, q in zip(psi, image, strict=True))
        residual = max(abs(q - lam * p).upper() for p, q in zip(psi, image, strict=True))
    else:
        transformed = fft(psi.conjugate() if fibonacci else psi)
        scale = (1 if fibonacci else cmath.exp(1j * cmath.pi * ((d - 1) % 24) / 12)) / np.sqrt(d)
        image = scale * np.exp(1j * np.pi * np.array(exponents) / d) * transformed[(-np.arange(d)) % d]
        lam = np.vdot(psi, image)
        residual = float(np.max(np.abs(image - lam * psi)))
    return residual, lam


def main():
    precise = "--precise" in sys.argv
    ctx.prec = 640  # About 192 decimal digits, with guard digits beyond the old 160.
    ctx.threads = 1
    files = [Path(arg) for arg in sys.argv[1:] if arg != "--precise"]
    if not files:
        files = [p for folder in ("fiducials", "fibonacci")
                 for p in sorted((Path(__file__).parent / folder).glob("fiducial_*.txt"))
                 if re.fullmatch(r"fiducial_\d+\.txt", p.name)]
    failed = False
    for path in files:
        rows = [re.findall(r"[+-]\.\d+e-?\d+", line)
                for line in path.read_text().splitlines() if line.strip()]
        parts = [(re_.replace(".", "0.", 1), im.replace(".", "0.", 1)) for re_, im in rows]
        psi = ([acb(re_, im) for re_, im in parts] if precise else
               np.array([complex(float(re_), float(im)) for re_, im in parts]))
        d = len(psi)
        fft = fourier(d, precise)
        worst = check(psi, precise, fft)
        fibonacci = path.resolve().parent.name == "fibonacci"
        residual, lam = symmetry(psi, fibonacci, precise, fft)
        named = re.fullmatch(r"fiducial_(\d+)\.txt", path.name)
        ok = d > 1 and (not named or int(named[1]) == d) and psi[0].imag == 0 and psi[0].real > 0
        ok = ok and max(worst, residual) < (arb("1e-147") if precise else 1e-12)
        failed |= not ok
        k = round(cmath.phase(complex(lam)) / (2 * cmath.pi / 3)) % 3
        if fibonacci:
            print(f"{path.name}: d = {d}, max deviation {float(worst):.1e}, Fibonacci F_f residual {float(residual):.1e} {'PASS' if ok else 'FAIL'}")
        else:
            print(f"{path.name}: d = {d}, max deviation {float(worst):.1e}, Zauner residual {float(residual):.1e}, eigenvalue ω^{k} {'PASS' if ok else 'FAIL'}")
    sys.exit(failed)


if __name__ == "__main__":
    main()
