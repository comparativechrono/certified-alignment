#!/usr/bin/env python3
"""Fetch exact official commits and build the native external comparison adapter."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "edlib": {
        "url": "https://github.com/Martinsos/edlib.git",
        "commit": "0ddc23ea06abd95ab6227442acb4c27bd7607b68",
        "description": "1.3.9.post1-8-g0ddc23e",
        "license": "MIT",
    },
    "WFA2-lib": {
        "url": "https://github.com/smarco/WFA2-lib.git",
        "commit": "3282f834713cb9767884949143fd0c803acae697",
        "description": "v2.3.6-21-g3282f83",
        "license": "MIT",
    },
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deps", type=Path, default=ROOT / "external_deps")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    deps = args.deps.resolve()
    deps.mkdir(parents=True, exist_ok=True)
    build = ROOT / "build"
    build.mkdir(exist_ok=True)
    commands = []

    def run(command, **kwargs):
        command = list(map(str, command))
        commands.append(command)
        subprocess.run(command, check=True, **kwargs)

    metadata = {}
    for name, spec in DEPENDENCIES.items():
        directory = deps / name
        if not directory.exists():
            run(["git", "clone", spec["url"], directory])
            run(["git", "-C", directory, "checkout", "--detach", spec["commit"]])
        actual = subprocess.check_output(
            ["git", "-C", str(directory), "rev-parse", "HEAD"], text=True).strip()
        if actual != spec["commit"]:
            raise RuntimeError(f"{name}: existing checkout is not the pinned commit")
        dirty = subprocess.check_output(
            ["git", "-C", str(directory), "diff", "--name-only", "HEAD"], text=True).strip()
        if dirty:
            raise RuntimeError(f"{name}: dependency checkout has modified tracked files")
        metadata[name] = dict(spec, license_sha256=hashlib.sha256(
            (directory / "LICENSE").read_bytes()).hexdigest())

    cxx = os.environ.get("CXX", "g++")
    cc = os.environ.get("CC", "gcc")
    flags = ["-O3", "-DNDEBUG", "-std=c++17"]
    edlib = deps / "edlib"
    wfa = deps / "WFA2-lib"
    # Edlib is one C++ translation unit; direct compilation avoids requiring CMake/Meson.
    run([cxx, *flags, "-I" + str(edlib / "edlib/include"), "-c",
         edlib / "edlib/src/edlib.cpp", "-o", build / "external_edlib.o"])
    # Rename only the old CLI entry point at compilation; engine/predictor bytes stay frozen.
    run([cxx, *flags, "-I" + str(ROOT / "include"), "-Dmain=certalign_cli_unused_main",
         "-c", ROOT / "src/aligner.cpp", "-o", build / "external_prototype.o"])
    # Command-line CC_FLAGS overrides WFA2's default -march=native for a common ISA baseline.
    # Force rebuild so an older build with different flags cannot be silently reused.
    run(["make", "-C", wfa, "-B", f"-j{args.jobs}", "BUILD_EXAMPLES=0", "BUILD_TOOLS=0",
         "BUILD_WFA_PARALLEL=0", "CC=" + cc, "CXX=" + cxx,
         "CC_FLAGS=-O3 -DNDEBUG -fPIE"])
    binary = build / "external_adapter"
    run([cxx, *flags, "-Wall", "-Wextra", "-I" + str(edlib / "edlib/include"),
         "-I" + str(wfa), "-I" + str(ROOT / "include"), ROOT / "benchmarks/external_adapter.cpp",
         build / "external_edlib.o", build / "external_prototype.o",
         wfa / "lib/libwfa.a", "-lm", "-o", binary])
    manifest = {
        "dependencies": metadata,
        "commands": commands,
        "cxx_version": subprocess.check_output([cxx, "--version"], text=True).splitlines()[0],
        "cc_version": subprocess.check_output([cc, "--version"], text=True).splitlines()[0],
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
        "adapter_source_sha256": hashlib.sha256(
            (ROOT / "benchmarks/external_adapter.cpp").read_bytes()).hexdigest(),
        "prototype_source_sha256": {str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted([ROOT / "src/aligner.cpp", *(ROOT / "include").glob("*.hpp")])},
        "common_optimization": "-O3 -DNDEBUG; no -march=native; single thread",
    }
    (build / "external_build.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"binary": str(binary), "manifest": str(build / "external_build.json")}))


if __name__ == "__main__":
    main()
