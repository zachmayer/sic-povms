# /// script
# requires-python = "==3.13.*"
# dependencies = ["numpy==2.5.4", "scipy==1.18.1"]
# [tool.uv]
# exclude-newer = "2026-10-19T00:00:00Z"
# ///
"""Recover a SIC seed from stark_large.jl's conjugated orbit of overlaps.

Unlike stark.jl's dense projector, FFT rows apply the partial projector to a
vector. Power iteration seeks a dominant positive eigenvalue; alternating
projections complete the missing overlaps. Rejecting other eigenpairs is a
search heuristic, not a proof that the corresponding shift cannot give a SIC.

Read metadata.json, ComplexF64 units.bin in Julia order, and Int64 label pairs.
The owner/sign maps use 5*d*d bytes of RAM, then read-only files. Each
numerical stage resumes from atomic files under --work. Representative choices
are in the Julia header.
"""

import argparse
import cmath
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time

for variable in ("OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np

from search import write_seed


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
    """Canonical fit ranking; pin rationale is in representatives.toml."""
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
                              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()+Path(__file__).with_name("search.py").read_bytes()).hexdigest()}
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
        """Map each standard displacement to a 1-based orbit index and its ±1 sign.

        Unknown entries have owner=sign=0; signs account for reducing labels mod d.
        The maps are built in RAM, which takes minutes, and written once. Scattered
        writes through a writable memmap made macOS write back many times the maps'
        size. An interrupted build starts again from zero.
        """
        checkpoint = self.work/"owner_progress.json"
        paths = (self.work/"owner.i32", self.work/"sign.i8")
        if not (checkpoint.exists() and json.loads(checkpoint.read_text())["complete"]):
            flat_owner = np.zeros(self.d*self.d, dtype=np.int32)
            flat_sign = np.zeros(self.d*self.d, dtype=np.int8)
            transform = np.asarray(self.metadata["standard_map"], dtype=np.int64).reshape(2, 2)
            matrices = [(transform @ np.asarray(s, dtype=np.int64).reshape(2, 2)) % self.dd
                        for s in self.metadata["stabilizers"]]
            started = time.perf_counter()
            for start in range(0, self.size, 262144):
                stop = min(start+262144, self.size)
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
                emit("owner_progress", groups=stop, total=self.size, seconds=time.perf_counter()-started)
            # Flat indices x+d*y are the column-major bytes of the (d, d) maps.
            for path, values in zip(paths, (flat_owner, flat_sign), strict=True):
                temporary = path.with_suffix(path.suffix+".part")
                with temporary.open("wb") as handle:
                    values.tofile(handle)
                    handle.flush()
                    os.fsync(handle.fileno())
                temporary.replace(path)
            save_json(checkpoint, {"next_group": self.size, "complete": True})
        owner, signs = (np.memmap(path, dtype=dtype, mode="r", shape=(self.d, self.d), order="F")
                        for path, dtype in zip(paths, ("int32", "int8"), strict=True))
        return owner, signs

    def determinant_shifts(self) -> np.ndarray:
        """First shift per det(M): labels give [a-q_b*b,b] for M=a*I+b*L.
        Shifts with det=1 are Clifford unitaries commuting with U_Z.
        """
        output = self.work/"determinant_shifts.npy"
        if output.exists():
            return np.load(output)
        first = np.full(self.dd, self.size, dtype=np.int64)
        qb = int(self.metadata["q_b"])
        for start in range(0, self.size, 262144):
            stop = min(start+262144, self.size)
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
        """Conjugate units at owner-shift in the product of cyclic quotient axes.

        Owners are 1-based; unknown owners may wrap because their signs are zero.
        """
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
        for start in range(state["next_shift"], len(self.shifts), 128):
            stop = min(start+128, len(self.shifts))
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
        for start in range(0, self.d, 16):
            yield np.arange(start, min(start+16, self.d), dtype=np.int64)

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
                if self.metadata["preferred_shift"] >= 0:
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
                residual = float("inf")
                eigenvalue = 0.0
                for iteration in range(state["iteration"], 100):
                    started = time.perf_counter()
                    image = self.partial_action(vector, shift)
                    eigenvalue = float(np.vdot(vector, image).real)
                    residual = float(np.linalg.norm(image-eigenvalue*vector))
                    emit("eigen_iteration", trial=trial+1, iteration=iteration,
                         eigenvalue=eigenvalue, residual=residual, seconds=time.perf_counter()-started)
                    if residual < 1e-10:
                        break
                    vector = normalized(image)
                    state["iteration"] = iteration+1
                    self.checkpoint(state, vector)
                emit("eigen_complete", trial=trial+1, eigenvalue=eigenvalue, residual=residual, threshold=threshold)
                if residual >= 1e-10 or eigenvalue <= threshold:
                    state = {"trial": trial+1, "stage": "new"}
                    save_json(progress, state)
                    continue
                mode = int(round(cmath.phase(np.vdot(vector, zauner(vector)))/(2*np.pi/3))) % 3
                vector = symmetry_project(vector, mode)
                state.update(stage="completion", iteration=0, eigenvalue=eigenvalue,
                             eigen_residual=residual, zauner_mode=mode)
                self.checkpoint(state, vector)
            for iteration in range(state["iteration"], 250):
                started = time.perf_counter()
                error, fixed, free = self.completion_action(vector, shift)
                emit("completion_iteration", trial=trial+1, iteration=iteration, error=error,
                     seconds=time.perf_counter()-started)
                if error < 1e-13:
                    state.update(stage="complete", max_overlap_error=error)
                    self.checkpoint(state, vector)
                    return vector
                vector = symmetry_project(fixed if error > 1e-9 else free, state["zauner_mode"])
                state["iteration"] = iteration+1
                self.checkpoint(state, vector)
            state = {"trial": trial+1, "stage": "new"}
            save_json(progress, state)
        raise ArithmeticError("No tested determinant representative produced a numerical SIC")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("export", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started, cpu_started = time.perf_counter(), time.process_time()
    reconstruction = Reconstruction(args)
    vector = reconstruction.run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix(args.output.suffix+".part")
    write_seed(temporary, vector)
    temporary.replace(args.output)
    moved = zauner(vector)
    lam = np.vdot(vector, moved)
    result = {"d": len(vector),
              "output": str(args.output), "seconds": time.perf_counter()-started,
              "cpu_seconds": time.process_time()-cpu_started, "peak_rss_bytes": peak_rss(),
              "norm_error": abs(float(np.vdot(vector, vector).real)-1),
              "zauner_residual": float(np.max(np.abs(moved-lam*vector))),
              "configuration": reconstruction.configuration}
    save_json(args.work/"result.json", result)
    emit("complete", **result)


if __name__ == "__main__":
    main()
