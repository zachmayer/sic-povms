# /// script
# requires-python = ">=3.13"
# dependencies = ["modal==1.6.2"]
# ///
"""Fresh Linux rebuilds, one container per published file, one remote function.

modal run --detach modal_rebuild.py
modal run --detach modal_rebuild.py --files fiducials/fiducial_5779.txt
modal volume get sic-large-rebuild /RUN_ID ./modal_results

The entrypoint prints RUN_ID and call IDs, then exits. Results, logs, stage times
and SHA-256 comparisons live on the Volume, independently of the client. A long
file gets sequential 22-hour calls, resuming only its own numerical checkpoints.
SIC_REBUILD_MEMORY_MIB defaults to 65536 MiB per file.
SIC_REBUILD_BLAS_THREADS (default 8) sets Julia's BLAS threads; the Python stages always pin one.
This uploads the local repository snapshot; no GitHub checkout or saved seed is
used. Only the selected output is removed from each isolated fresh working copy.
Published large-d files use representatives.toml's pins from an unpublished
prototype run so its Mac runs and Linux rebuilds compare byte for byte. These
rebuilds exercise pinned choices; test_large.py tests the canonical rule at small d.
"""

import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import time
import uuid

import modal


ROOT = Path(__file__).resolve().parent
DEFAULT_FILES = (
    "fiducials/fiducial_2208.txt",
    "fibonacci/fiducial_2208.txt",
    "fiducials/fiducial_5779.txt",
    "fiducials/fiducial_15128.txt",
)
JULIA_URL = "https://julialang-s3.julialang.org/bin/linux/x64/1.13/julia-1.13.2-linux-x86_64.tar.gz"
VOLUME_NAME = "sic-large-rebuild"
COPY_IGNORE = (".git", "large_d_work", "modal_results", "__pycache__", "*.checkpoint.json*",
               "*.jsonl", "fiducial_*_stark.txt", "fiducial_*_seed*.txt", "fiducial_48g.txt")


def commands_for(filename: str, workers: int = 8) -> list[list[str]]:
    """The repository's construction commands, with explicit output ownership."""
    match = re.fullmatch(r"(fiducials|fibonacci)/fiducial_([1-9][0-9]*)\.txt", filename)
    if match is None or not 1 <= workers <= 8:
        raise ValueError("Expected a published fiducial path and 1..8 workers")
    folder, dimension = match.groups()
    if folder == "fibonacci":
        if dimension != "2208":
            raise ValueError("Only the published Fibonacci construction is supported")
        construction = ["uv", "run", "fibonacci_2208.py", "--fetch-input"]
        seed = f"fiducial_{dimension}_fibonacci_seed.txt"
    else:
        construction = ["julia", "--project=.", "--threads=1", "--gcthreads=1", "stark_large.jl", dimension]
        seed = f"fiducial_{dimension}_stark.txt"
    return [construction, ["uv", "run", "refine_large.py", seed, "--output", filename,
                           "--workers", str(workers)]]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_fingerprint(repo: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(repo.iterdir()):
        if path.suffix in (".py", ".jl", ".toml"):
            digest.update(path.name.encode("ascii") + b"\0" + path.read_bytes())
    return digest.hexdigest()


def save_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(".pending")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="ascii")
    temporary.replace(path)


def run_command(command: list[str], cwd: Path, log: Path, environment: dict[str, str],
                deadline: float, commit) -> int | None:
    """Stream logs and commit checkpoints; None means a clean slice interruption."""
    with log.open("ab") as stream:
        stream.write(("\nCOMMAND " + json.dumps(command) + "\n").encode("ascii"))
        stream.flush()
        child = subprocess.Popen(command, cwd=cwd, env=environment, stdout=stream,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        next_commit = time.monotonic() + 60
        while child.poll() is None:
            if time.monotonic() >= deadline:
                os.killpg(child.pid, signal.SIGTERM)
                try:
                    child.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(child.pid, signal.SIGKILL)
                    child.wait()
                commit()
                return None
            if time.monotonic() >= next_commit:
                commit()
                next_commit = time.monotonic() + 60
            time.sleep(0.1)
    commit()
    return child.returncode


def rebuild_file(repo: Path, destination: Path, filename: str, expected: str,
                 environment: dict[str, str], commit, *, workers: int = 8,
                 slice_seconds: float = 22 * 3600, timed: bool = True) -> dict:
    """Service-independent implementation, exercised locally on small dimensions."""
    commands = commands_for(filename, workers)
    destination.mkdir(parents=True, exist_ok=True)
    state_path = destination / "result.json"
    work = destination / "repo"
    source = source_fingerprint(repo)
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if (state["file"], state["committed_sha256"], state["commands"]) != (filename, expected, commands):
            raise ValueError("Rebuild identity changed")
        if state["source_sha256"] != source:
            raise ValueError("Rebuild source changed; start a new run")
    else:
        shutil.copytree(repo, work, ignore=shutil.ignore_patterns(*COPY_IGNORE))
        (work / filename).unlink(missing_ok=True)
        state = {"file": filename, "committed_sha256": expected, "commands": commands,
                 "source_sha256": source, "platform": platform.platform(),
                 "blas_threads": environment.get("SIC_BLAS_THREADS", "unspecified"),
                 "next_command": 0, "attempts": 0, "stages": [], "status": "RUNNING"}
    state["attempts"] += 1
    save_json(state_path, state)
    commit()
    deadline = time.monotonic() + slice_seconds
    for index in range(state["next_command"], len(commands)):
        command = commands[index]
        if timed:
            command = ["/usr/bin/time", "-v", *command]
        started = time.monotonic()
        result = run_command(command, work, destination / f"stage_{index}.log", environment, deadline, commit)
        state["stages"].append({"command": index, "seconds": time.monotonic() - started,
                                "exit_code": result})
        if result is None:
            state["status"] = "CONTINUE"
            save_json(state_path, state)
            commit()
            return state
        if result:
            state["status"] = "COMMAND_FAILED"
            save_json(state_path, state)
            commit()
            return state
        state["next_command"] = index + 1
        save_json(state_path, state)
        commit()
    rebuilt = work / filename
    state["rebuilt_sha256"] = sha256(rebuilt)
    state["status"] = "BYTE_IDENTICAL" if state["rebuilt_sha256"] == expected else "DIFFERENT"
    shutil.copyfile(rebuilt, destination / rebuilt.name)
    save_json(state_path, state)
    commit()
    return state


image = (
    modal.Image.debian_slim(python_version="3.13")
    .apt_install("curl", "time", "ca-certificates")
    .run_commands(f"curl -fsSL {JULIA_URL} | tar -xz -C /opt && ln -s /opt/julia-1.13.2/bin/julia /usr/local/bin/julia")
    .run_commands("curl -LsSf https://astral.sh/uv/install.sh | env UV_INSTALL_DIR=/usr/local/bin sh")
    .add_local_dir(ROOT, "/source", copy=True, ignore=list(COPY_IGNORE))
    .run_commands("cd /source && julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'")
)
app = modal.App("sic-large-rebuild", image=image)
volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)


@app.function(cpu=8, memory=int(os.environ.get("SIC_REBUILD_MEMORY_MIB", "65536")),
              timeout=24 * 3600, single_use_containers=True, volumes={"/results": volume})
def rebuild(run_id: str, filename: str, expected: str, blas_threads: int) -> dict:
    volume.reload()
    destination = Path("/results") / run_id / filename.replace("/", "_").removesuffix(".txt")
    environment = os.environ.copy()
    environment.update(SIC_BLAS_THREADS=str(blas_threads), OPENBLAS_NUM_THREADS=str(blas_threads),
                       OMP_NUM_THREADS=str(blas_threads), JULIA_NUM_THREADS="1")
    state = rebuild_file(Path("/source"), destination, filename, expected, environment, volume.commit)
    print(json.dumps(state), flush=True)
    if state["status"] == "CONTINUE":
        if state["attempts"] >= 8:
            state["status"] = "SLICE_LIMIT"
            save_json(destination / "result.json", state)
            volume.commit()
        else:
            rebuild.spawn(run_id, filename, expected, blas_threads)
    return state


@app.local_entrypoint()
def main(files: str = ",".join(DEFAULT_FILES)) -> None:
    names = files.split(",")
    if len(set(names)) != len(names):
        raise ValueError("Each published file must have exactly one writer")
    for name in names:
        commands_for(name)
        if not (ROOT / name).is_file():
            raise ValueError(f"No committed file to compare: {name}")
    run_id = uuid.uuid4().hex
    blas_threads = int(os.environ.get("SIC_REBUILD_BLAS_THREADS", "8"))
    if not 1 <= blas_threads <= 8:
        raise ValueError("SIC_REBUILD_BLAS_THREADS must be 1..8")
    print(f"RUN_ID={run_id}", flush=True)
    for name in names:
        call = rebuild.spawn(run_id, name, sha256(ROOT / name), blas_threads)
        print(f"{name}: {call.object_id}", flush=True)
    print(f"modal volume get {VOLUME_NAME} /{run_id} ./modal_results/{run_id}", flush=True)
