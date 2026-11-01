# /// script
# requires-python = ">=3.13"
# dependencies = ["modal==1.6.2"]
# ///
"""Offline small-d acceptance, without any cloud calls (at most two cores).

uv run test_large.py --work /absolute/path/outside/repo

Requires Julia on PATH and the repository's Julia/Python packages installed.
Runs fresh d=19 (odd, q=1) and d=48 (even, q=3), verifies their saved digits
independently, then exercises the Modal worker's actual logic with interruption,
continuation, fresh construction and committed-hash comparison. No large-d work.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from modal_rebuild import COPY_IGNORE, commands_for, rebuild_file, run_command


ROOT = Path(__file__).resolve().parent


def checked(command: list[str], cwd: Path, environment: dict[str, str], log: Path) -> None:
    with log.open("wb") as stream:
        subprocess.run(command, cwd=cwd, env=environment, stdout=stream,
                       stderr=subprocess.STDOUT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    if work == ROOT or ROOT in work.parents:
        parser.error("Keep acceptance artifacts outside the repository")
    work.mkdir(parents=True, exist_ok=False)
    environment = os.environ.copy()
    environment.update(SIC_BLAS_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
                       JULIA_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", MKL_NUM_THREADS="1")
    fixture = work / "fixture"
    shutil.copytree(ROOT, fixture, ignore=shutil.ignore_patterns(*COPY_IGNORE))
    results = {}
    for d in (19, 48):
        filename = f"fiducials/fiducial_{d}.txt"
        for index, command in enumerate(commands_for(filename, workers=1)):
            checked(command, fixture, environment, work / f"d{d}_baseline_{index}.log")
        checked(["uv", "run", "verify.py", filename, "--precise"], fixture,
                environment, work / f"d{d}_precise.log")
        expected = hashlib.sha256((fixture / filename).read_bytes()).hexdigest()
        destination = work / f"rebuild_{d}"
        # Interrupt the actual command supervisor before construction completes.
        interrupted = rebuild_file(fixture, destination, filename, expected, environment,
                                   lambda: None, workers=2, slice_seconds=0.01, timed=False)
        assert interrupted["status"] == "CONTINUE", interrupted
        state = rebuild_file(fixture, destination, filename, expected, environment,
                             lambda: None, workers=2, timed=False)
        assert state["status"] == "BYTE_IDENTICAL", state
        # Repeated delivery reads only completed state and hashes the output.
        repeated = rebuild_file(fixture, destination, filename, expected, environment,
                                lambda: None, workers=2, timed=False)
        assert repeated["status"] == "BYTE_IDENTICAL", repeated
        results[str(d)] = state
    # A wrong expected hash must produce DIFFERENT, never successful acceptance.
    destination = work / "rebuild_19"
    path = destination / "result.json"
    state = json.loads(path.read_text())
    state["committed_sha256"] = "0" * 64
    path.write_text(json.dumps(state))
    different = rebuild_file(fixture, destination, "fiducials/fiducial_19.txt", "0" * 64,
                             environment, lambda: None, workers=2, timed=False)
    assert different["status"] == "DIFFERENT", different
    failure = run_command(["uv", "run", "--no-project", "python", "-c", "raise SystemExit(7)"],
                          fixture, work / "failure.log", environment, float("inf"), lambda: None)
    assert failure == 7, failure
    for invalid in ("../fiducial_19.txt", "/fiducials/fiducial_19.txt", "fibonacci/fiducial_48.txt"):
        try:
            commands_for(invalid)
        except ValueError:
            continue
        raise AssertionError(f"Accepted invalid rebuild target: {invalid}")
    (work / "acceptance.json").write_text(json.dumps({"status": "PASS", "results": results,
                                                     "negative_controls": "PASS"}, indent=2) + "\n")
    print("PASS: fresh d=19 and q=3 d=48, precise verification, two-worker byte equality, restart and negative controls")


if __name__ == "__main__":
    main()
