# /// script
# requires-python = "==3.13.*"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1", "python-flint==0.9.0"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Stream partial-projector SIC reconstruction without dense matrices.

The input directory contains metadata.json, complex128 units.bin in Julia/F
order, and int64 labels.bin as N consecutive (p1,p2) pairs. Metadata explicitly
provides d, axes, q_b, stabilizers (flat row-major 2x2), and standard_map.
Owner/sign maps cost 5*d*d bytes on disk and are memory mapped. No complex d*d
array is allocated. Every completed owner chunk, fit batch, eigen iteration,
and alternating-projection iteration can resume after interruption.

The default power iteration uses the largest-magnitude eigenpair. A genuine
partial rank-one projector is expected to have a dominant positive eigenvalue;
negative candidates are rejected. This expectation is a numerical search
heuristic, not a proof that all unsuccessful shifts can be excluded.
"""

import argparse
import cmath
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import os
for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

from search import fiducial_text
from fractions import Fraction

def write_fiducial(path: Path, psi: np.ndarray) -> None:
    psi = psi / np.linalg.norm(psi)
    psi *= np.conj(psi[0]) / abs(psi[0])
    path.write_bytes(fiducial_text([(Fraction(float(z.real)), Fraction(0) if j == 0 else Fraction(float(z.imag)))
                                  for j, z in enumerate(psi)], 17).encode("ascii"))


def emit(event: str, **fields: object) -> None:
    print(json.dumps({"event": event, **fields}), flush=True)


def save_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix+".part")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n")
    temporary.replace(path)


def save_array(path: Path, value: np.ndarray) -> None:
    temporary = path.with_suffix(path.suffix+".part")
    with temporary.open("wb") as handle:
        np.save(handle, value)
    temporary.replace(path)


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8*1024*1024), b""):
            result.update(block)
    return result.hexdigest()


def peak_rss() -> int:
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if sys.platform == "darwin" else value*1024)


def chirp(d: int, shifts: np.ndarray) -> np.ndarray:
    exponents = ((d+1)*shifts[:, None]*np.arange(d, dtype=np.int64)) % (2*d)
    return np.exp(1j*np.pi*exponents/d)


def normalized(vector: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(vector)
    if not np.isfinite(norm) or norm == 0:
        raise ArithmeticError("Nonfinite or zero reconstruction vector")
    return vector/norm


def shift_order(fits: np.ndarray, shifts: np.ndarray, preferred: int) -> np.ndarray:
    """Canonical fit ranking, or the pinned existing Clifford representative."""
    if preferred < 0:
        return np.argsort(np.round(fits, 10), kind="stable")
    order = np.flatnonzero(shifts == preferred)
    if len(order) != 1:
        raise ValueError("Pinned representative is not a determinant-class shift")
    return order


def zauner(vector: np.ndarray) -> np.ndarray:
    d = len(vector)
    r = np.arange(d, dtype=np.int64)
    phases = np.exp(1j*np.pi*(((d+1)*r*r) % (2*d))/d)
    return cmath.exp(1j*np.pi*((d-1) % 24)/12)*np.sqrt(d)*phases*np.fft.ifft(vector)


def symmetry_project(vector: np.ndarray, mode: int) -> np.ndarray:
    lam = cmath.exp(2j*np.pi*mode/3)
    moved = zauner(vector)
    return normalized((vector+lam.conjugate()*moved+lam.conjugate()**2*zauner(moved))/3)


def add_projector_rows(output: np.ndarray, vector: np.ndarray, shifts: np.ndarray,
                       rows: np.ndarray) -> None:
    """Add (Braw+Braw*)/2; caller adds the identity/d once."""
    d = len(vector)
    indices = (shifts[:, None]+np.arange(d, dtype=np.int64)) % d
    output += np.sum(np.conj(rows)*vector[indices], axis=0)/2
    for shift, row in zip(shifts, rows, strict=True):
        output += np.roll(row*vector, int(shift))/2


class Reconstruction:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.work = args.work
        self.work.mkdir(parents=True, exist_ok=True)
        metadata_path = args.export/"metadata.json"
        self.metadata = json.loads(metadata_path.read_text())
        self.d = int(self.metadata["d"])
        self.dd = self.d*(2 if self.d % 2 == 0 else 1)
        self.axes = tuple(int(x) for x in self.metadata["axes"])
        self.size = int(np.prod(self.axes, dtype=np.int64))
        if self.size >= np.iinfo(np.int32).max:
            raise ValueError("Group exceeds int32 owner capacity")
        for name, size in (("units.bin", 16*self.size), ("labels.bin", 16*self.size)):
            if (args.export/name).stat().st_size != size:
                raise ValueError(f"Incorrect {name} size: expected {size}")
        self.configuration = {"metadata_sha256": digest(metadata_path),
                              "units_sha256": digest(args.export/"units.bin"),
                              "labels_sha256": digest(args.export/"labels.bin"),
                              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name("search.py").read_bytes()).hexdigest(),
                              "eigen_method": args.eigen_method,
                              "eigen_tolerance": args.eigen_tolerance,
                              "eigen_steps": args.eigen_steps,
                              "completion_steps": args.completion_steps,
                              "tolerance": args.tolerance}
        manifest = self.work/"configuration.json"
        if manifest.exists():
            if json.loads(manifest.read_text()) != self.configuration:
                raise ValueError("Reconstruction checkpoint input/code/settings mismatch")
        else:
            save_json(manifest, self.configuration)
        self.units = np.memmap(args.export/"units.bin", dtype="<c16", mode="r", shape=(self.size,))
        self.labels = np.memmap(args.export/"labels.bin", dtype="<i8", mode="r", shape=(self.size, 2))
        self.owner, self.sign = self.build_owner()
        self.known_count = int(np.count_nonzero(self.owner))+int(self.owner[0, 0] == 0)
        self.shifts = self.determinant_shifts()

    def build_owner(self) -> tuple[np.memmap, np.memmap]:
        checkpoint = self.work/"owner_progress.json"
        state = json.loads(checkpoint.read_text()) if checkpoint.exists() else {"next_group": 0}
        mode = "r+" if checkpoint.exists() else "w+"
        owner = np.memmap(self.work/"owner.i32", dtype="int32", mode=mode, shape=(self.d, self.d), order="F")
        signs = np.memmap(self.work/"sign.i8", dtype="int8", mode=mode, shape=(self.d, self.d), order="F")
        flat_owner, flat_sign = owner.reshape(-1, order="F"), signs.reshape(-1, order="F")
        transform = np.asarray(self.metadata["standard_map"], dtype=np.int64).reshape(2, 2)
        matrices = [(transform @ np.asarray(s, dtype=np.int64).reshape(2, 2)) % self.dd
                    for s in self.metadata["stabilizers"]]
        started = time.perf_counter()
        for start in range(state["next_group"], self.size, self.args.owner_chunk):
            stop = min(start+self.args.owner_chunk, self.size)
            labels = self.labels[start:stop]
            owners = np.arange(start+1, stop+1, dtype=np.int32)
            for matrix in matrices:
                a = (matrix[0, 0]*labels[:, 0]+matrix[0, 1]*labels[:, 1]) % self.dd
                b = (matrix[1, 0]*labels[:, 0]+matrix[1, 1]*labels[:, 1]) % self.dd
                x, y = a % self.d, b % self.d
                indices = x+self.d*y
                # Preserve the original Julia g-then-stabilizer overwrite order.
                np.maximum.at(flat_owner, indices, owners)
                winners = owners == flat_owner[indices]
                parity = ((a*b-x*y)//self.d) % 2
                flat_sign[indices[winners]] = 1-2*parity[winners]
            owner.flush()
            signs.flush()
            save_json(checkpoint, {"next_group": stop, "complete": stop == self.size})
            emit("owner_progress", groups=stop, total=self.size, seconds=time.perf_counter()-started)
        return owner, signs

    def determinant_shifts(self) -> np.ndarray:
        output = self.work/"determinant_shifts.npy"
        if output.exists():
            return np.load(output)
        first = np.full(self.dd, self.size, dtype=np.int64)
        qb = int(self.metadata["q_b"])
        for start in range(0, self.size, self.args.owner_chunk):
            stop = min(start+self.args.owner_chunk, self.size)
            labels = self.labels[start:stop]
            a = (labels[:, 0]+qb*labels[:, 1]) % self.dd
            b = labels[:, 1]
            determinant = (a*a-qb*a*b+b*b) % self.dd
            np.minimum.at(first, determinant, np.arange(start, stop, dtype=np.int64))
        result = np.sort(first[first < self.size])
        save_array(output, result)
        emit("determinant_shifts", count=len(result), known=self.known_count, total=self.d*self.d)
        return result

    def shifted_units(self, owners: np.ndarray, shift: int | np.ndarray) -> np.ndarray:
        flat = owners.astype(np.int64)-1
        result = np.zeros(np.broadcast_shapes(flat.shape, np.shape(shift)), dtype=np.int64)
        stride = 1
        for axis in self.axes:
            result += (((flat//stride) % axis-(shift//stride) % axis) % axis)*stride
            stride *= axis
        return np.conj(self.units[result])

    def known_rows(self, shifts: np.ndarray, group_shift: int) -> np.ndarray:
        owners = self.owner[shifts, :]
        return self.sign[shifts, :]*self.shifted_units(owners, group_shift)

    def fit_shifts(self) -> np.ndarray:
        path = self.work/"shift_fits.npy"
        progress = self.work/"fit_progress.json"
        state = json.loads(progress.read_text()) if progress.exists() else {"next_shift": 0}
        fits = np.load(path) if path.exists() else np.full(len(self.shifts), np.nan)
        tau = chirp(self.d, np.array([0, 1]))/(self.d+1)
        started = time.perf_counter()
        for start in range(state["next_shift"], len(self.shifts), self.args.shift_batch):
            stop = min(start+self.args.shift_batch, len(self.shifts))
            shifts = self.shifts[start:stop, None]
            diagonal = np.fft.ifft(self.sign[0, :]*self.shifted_units(self.owner[0, :], shifts)*tau[0], axis=1)
            diagonal = 1/self.d+diagonal.real
            off_diagonal = np.fft.ifft(self.sign[1, :]*self.shifted_units(self.owner[1, :], shifts)*tau[1], axis=1)
            fits[start:stop] = self.d*np.sum(np.abs(np.abs(off_diagonal)**2-diagonal*np.roll(diagonal, -1, axis=1)), axis=1)
            save_array(path, fits)
            save_json(progress, {"next_shift": stop, "complete": stop == len(self.shifts)})
            emit("fit_progress", shifts=stop, total=len(self.shifts), seconds=time.perf_counter()-started)
        if not np.isfinite(fits).all():
            raise ArithmeticError("Incomplete or nonfinite shift fits")
        return fits

    def batches(self):
        for start in range(0, self.d, self.args.row_batch):
            yield np.arange(start, min(start+self.args.row_batch, self.d), dtype=np.int64)

    def partial_action(self, vector: np.ndarray, shift: int) -> np.ndarray:
        output = vector/self.d
        for a in self.batches():
            rows = np.fft.ifft(self.known_rows(a, shift)*chirp(self.d, a), axis=1)/(self.d+1)
            add_projector_rows(output, vector, a, rows)
        return output

    def completion_action(self, vector: np.ndarray, shift: int) -> tuple[float, np.ndarray, np.ndarray]:
        norm = np.vdot(vector, vector).real
        fixed, free = vector/self.d, vector/self.d
        maximum = 0.0
        for a in self.batches():
            indices = (a[:, None]+np.arange(self.d, dtype=np.int64)) % self.d
            tau = chirp(self.d, a)
            overlaps = np.conj((self.d+1)*self.d/norm*tau*np.fft.ifft(np.conj(vector[indices])*vector, axis=1))
            errors = np.abs(np.abs(overlaps)**2-(self.d+1))/(self.d+1)**2
            if a[0] == 0:
                errors[0, 0] = 0
            maximum = max(maximum, float(np.max(errors)))
            moduli = np.abs(overlaps)
            if np.any(moduli == 0) or not np.isfinite(moduli).all():
                raise ArithmeticError("Completion encountered zero/nonfinite overlap; phase projection undefined")
            target = np.sqrt(self.d+1)*overlaps/moduli
            rows = np.fft.ifft(target*tau, axis=1)/(self.d+1)
            add_projector_rows(free, vector, a, rows)
            known = self.owner[a, :] != 0
            if a[0] == 0:
                known[0, 0] = True
            rows = np.fft.ifft(np.where(known, self.known_rows(a, shift), target)*tau, axis=1)/(self.d+1)
            add_projector_rows(fixed, vector, a, rows)
        return maximum, fixed, free

    def checkpoint(self, state: dict, vector: np.ndarray) -> None:
        generation = int(state.get("generation", 0))+1
        filename = f"vector_{generation % 2}.npy"
        save_array(self.work/filename, vector)
        state.update(generation=generation, vector_file=filename)
        save_json(self.work/"search_progress.json", state)

    def run(self) -> np.ndarray:
        fits = self.fit_shifts()
        order = shift_order(fits, self.shifts, int(self.metadata["preferred_shift"]))
        progress = self.work/"search_progress.json"
        state = json.loads(progress.read_text()) if progress.exists() else {"trial": 0, "stage": "new"}
        if state["stage"] == "complete":
            return np.load(self.work/state["vector_file"])
        threshold = self.known_count/self.d**2/2
        for trial in range(state["trial"], len(order)):
            selected = int(order[trial])
            shift = int(self.shifts[selected])
            if state["trial"] != trial or state["stage"] == "new":
                if self.metadata["legacy_start"]:
                    rng = np.random.default_rng(221+trial)
                    vector = normalized(rng.standard_normal(self.d)+1j*rng.standard_normal(self.d))
                else:
                    vector = np.zeros(self.d, dtype=complex)
                    vector[0] = 1
                state = {"trial": trial, "selected": selected, "group_shift": shift,
                         "stage": "eigen", "iteration": 0, "fit": float(fits[selected])}
                self.checkpoint(state, vector)
            else:
                vector = np.load(self.work/state["vector_file"])
            emit("shift_trial", trial=trial+1, selected=selected, group_shift=shift, fit=float(fits[selected]))
            if state["stage"] == "eigen":
                if self.args.eigen_method == "eigsh":
                    operator = LinearOperator((self.d, self.d), matvec=lambda v: self.partial_action(v, shift), dtype=np.complex128)
                    values, vectors = eigsh(operator, k=1, which="LA", v0=vector,
                                            tol=self.args.eigen_tolerance, maxiter=self.args.eigen_steps)
                    vector = normalized(vectors[:, 0])
                    image = self.partial_action(vector, shift)
                    eigenvalue = float(values[0])
                    residual = float(np.linalg.norm(image-eigenvalue*vector))
                else:
                    residual = float("inf")
                    eigenvalue = 0.0
                    for iteration in range(state["iteration"], self.args.eigen_steps):
                        started = time.perf_counter()
                        image = self.partial_action(vector, shift)
                        eigenvalue = float(np.vdot(vector, image).real)
                        residual = float(np.linalg.norm(image-eigenvalue*vector))
                        emit("eigen_iteration", trial=trial+1, iteration=iteration,
                             eigenvalue=eigenvalue, residual=residual, seconds=time.perf_counter()-started)
                        if residual < self.args.eigen_tolerance:
                            break
                        vector = normalized(image)
                        state["iteration"] = iteration+1
                        self.checkpoint(state, vector)
                emit("eigen_complete", trial=trial+1, eigenvalue=eigenvalue, residual=residual, threshold=threshold)
                if residual >= self.args.eigen_tolerance or eigenvalue <= threshold:
                    state = {"trial": trial+1, "stage": "new"}
                    save_json(progress, state)
                    continue
                mode = int(round(cmath.phase(np.vdot(vector, zauner(vector)))/(2*np.pi/3))) % 3
                vector = symmetry_project(vector, mode)
                state.update(stage="completion", iteration=0, eigenvalue=eigenvalue,
                             eigen_residual=residual, zauner_mode=mode)
                self.checkpoint(state, vector)
            for iteration in range(state["iteration"], self.args.completion_steps):
                started = time.perf_counter()
                error, fixed, free = self.completion_action(vector, shift)
                emit("completion_iteration", trial=trial+1, iteration=iteration, error=error,
                     seconds=time.perf_counter()-started)
                if error < self.args.tolerance:
                    state.update(stage="complete", max_overlap_error=error)
                    self.checkpoint(state, vector)
                    return vector
                vector = symmetry_project(fixed if error > 1e-9 else free, state["zauner_mode"])
                state["iteration"] = iteration+1
                self.checkpoint(state, vector)
            state = {"trial": trial+1, "stage": "new"}
            save_json(progress, state)
        raise ArithmeticError("No tested determinant representative produced a numerical SIC")


def self_test() -> None:
    rng = np.random.default_rng(39604)
    maximum = 0.0
    for d in (3, 4, 8, 13):
        w = rng.standard_normal((d, d))+1j*rng.standard_normal((d, d))
        vector = normalized(rng.standard_normal(d)+1j*rng.standard_normal(d))
        shifts = np.arange(d, dtype=np.int64)
        rows = np.fft.ifft(w*chirp(d, shifts), axis=1)/(d+1)
        matrix = np.array([[rows[(j-k) % d, k]+(j == k)/d for k in range(d)] for j in range(d)])
        matrix = (matrix+matrix.conj().T)/2
        image = vector/d
        for start in range(0, d, 3):
            add_projector_rows(image, vector, shifts[start:start+3], rows[start:start+3])
        error = float(np.max(np.abs(image-matrix@vector)))
        if error > 2e-14:
            raise ArithmeticError(f"Explicit projector mismatch at d={d}: {error}")
        maximum = max(maximum, error)
        for mode in range(3):
            lam = cmath.exp(2j*np.pi*mode/3)
            raw = vector+lam.conjugate()*zauner(vector)+lam.conjugate()**2*zauner(zauner(vector))
            if np.linalg.norm(raw) < 1e-12:
                continue  # Some dimensions have a zero-dimensional Zauner eigenspace.
            projected = symmetry_project(vector, mode)
            residual = np.max(np.abs(zauner(projected)-cmath.exp(2j*np.pi*mode/3)*projected))
            if residual > 2e-13:
                raise ArithmeticError(f"Zauner projection mismatch at d={d}, mode={mode}")
    emit("self_test", verdict="PASS", dimensions=[3, 4, 8, 13], maximum_projector_error=maximum)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path, nargs="?")
    parser.add_argument("--work", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--owner-chunk", type=int, default=262144)
    parser.add_argument("--shift-batch", type=int, default=128)
    parser.add_argument("--row-batch", type=int, default=16)
    parser.add_argument("--eigen-method", choices=("power", "eigsh"), default="power")
    parser.add_argument("--eigen-steps", type=int, default=100)
    parser.add_argument("--eigen-tolerance", type=float, default=1e-10)
    parser.add_argument("--completion-steps", type=int, default=250)
    parser.add_argument("--tolerance", type=float, default=1e-13)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.export or not args.work or not args.output:
        raise ValueError("Reconstruction requires EXPORT --work DIR --output FILE")
    if min(args.owner_chunk, args.shift_batch, args.row_batch, args.eigen_steps, args.completion_steps) < 1:
        raise ValueError("Chunk sizes and iteration limits must be positive")
    started, cpu_started = time.perf_counter(), time.process_time()
    reconstruction = Reconstruction(args)
    vector = reconstruction.run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix+".part")
    write_fiducial(temporary, vector)
    temporary.replace(args.output)
    lam = np.vdot(vector, zauner(vector))
    result = {"candidate_only_requires_high_precision_refinement": True, "d": len(vector),
              "output": str(args.output), "seconds": time.perf_counter()-started,
              "cpu_seconds": time.process_time()-cpu_started, "peak_rss_bytes": peak_rss(),
              "norm_error": abs(float(np.vdot(vector, vector).real)-1),
              "zauner_residual": float(np.max(np.abs(zauner(vector)-lam*vector))),
              "owner_and_sign_file_bytes": 5*len(vector)**2,
              "configuration": reconstruction.configuration}
    save_json(args.work/"result.json", result)
    emit("complete", **result)


if __name__ == "__main__":
    main()
