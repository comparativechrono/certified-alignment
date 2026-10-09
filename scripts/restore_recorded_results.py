#!/usr/bin/env python3
"""Restore losslessly compressed historical outputs to their original paths."""
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
entries = json.loads((ROOT / "results/COMPRESSED_RESULTS.json").read_text())["files"]
for entry in entries:
    packed = (ROOT / entry["compressed_path"]).read_bytes()
    if hashlib.sha256(packed).hexdigest() != entry["compressed_sha256"]:
        raise SystemExit(f"compressed checksum mismatch: {entry['compressed_path']}")
    raw = gzip.decompress(packed)
    if len(raw) != entry["bytes"] or hashlib.sha256(raw).hexdigest() != entry["sha256"]:
        raise SystemExit(f"original checksum mismatch: {entry['path']}")
    target = ROOT / entry["path"]
    if target.exists() and target.read_bytes() != raw:
        raise SystemExit(f"refusing to overwrite changed output: {entry['path']}")
    target.write_bytes(raw)
print(f"Verified and restored {len(entries)} historical result files.")
