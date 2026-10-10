#!/usr/bin/env python3
"""Fetch the checksum-pinned historical Locityper panel and selected loci.

Only regular, explicitly selected archive members are materialized. Byte-range
downloads are checked against Content-Range and the entire official MD5; SHA256
and member hashes are recorded for independently checking downstream analyses.
The 785 MB archive and extracted FASTA/PAF data are intentionally not Git assets.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import time
import urllib.request

SOURCE = "https://zenodo.org/records/15879519/files/db.tar.gz?download=1"
SIZE = 785474256
MD5 = "f782f6e72426b7a818c7e1fd42f50277"
CHUNK = 64 * 1024 * 1024
MEMBERS = {"ref.bed", "discarded_haplotypes.txt", "haplotypes.fa.gz", "haplotypes.paf.gz"}
EVAL_COMMIT = "146cfd42e9179bafd63d1f41db7fc449940989de"
EVAL_HASHES = {
    "extra/gt_dist.py": "1099c43e4489672d5362395e8fb6cdaa8842b32e5d925b38d8fa565b30a5c26d",
    "extra/eval_accuracy.py": "20f18314e27416c0a1f98447d7854993c9e402223e0174ff450f57f98da44eb6",
    "extra/common.py": "67e4c151f8e7c55efd0ef662ae24353357347816e08c61e680276bba6b3927d8",
    "LICENSE": "f90d1546fe0e436830365c9196856450a72896ee25472fe450dc800f8b7ef070",
}
REFERENCE_EVIDENCE = {
    "assembly": "GRCh38",
    "paper": "https://doi.org/10.1038/s41588-025-02362-4",
    "paper_fulltext_xml": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12597825/fullTextXML",
    "paper_evidence": "Data availability specifies the HPRC GRCh38 graph raw VCF; Supplementary Table 1 is explicitly labelled GRCh38 coordinates.",
    "crosscheck": "All 15 selected archived ref.bed intervals exactly match COSIGT cmrgs_original.bed, snapshotted in data/panel_adequacy/cosigt/.",
    "coordinate_convention": "0-based half-open BED. Use these original intervals, not COSIGT refined intervals.",
}


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as src:
        for block in iter(lambda: src.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def fetch_range(directory, start, end, source=SOURCE, size=SIZE):
    path = directory / f"{start:012d}.part"
    expected_size = end - start + 1
    for attempt in range(8):
        existing = path.stat().st_size if path.exists() else 0
        if existing == expected_size:
            return path
        if existing > expected_size:
            raise ValueError(f"Oversized range part: {path}")
        current = start + existing
        req = urllib.request.Request(source, headers={"Range": f"bytes={current}-{end}"})
        try:
            with urllib.request.urlopen(req, timeout=90) as response:
                wanted = f"bytes {current}-{end}/{size}"
                if response.status != 206 or response.headers.get("Content-Range") != wanted:
                    raise ValueError(f"Server did not honor range: {response.status}, {response.headers}")
                with path.open("ab") as out:
                    while block := response.read(256 * 1024):
                        out.write(block)
                    out.flush()
                    os.fsync(out.fileno())
            if path.stat().st_size == expected_size:
                print(f"Downloaded range {start}-{end}", flush=True)
                return path
        except (OSError, TimeoutError) as exc:
            print(f"Retry range {start}-{end}: {exc}", flush=True)
        time.sleep(min(attempt + 1, 5))
    raise RuntimeError(f"Incomplete range after retries: {path}")


def download(root, workers, source=SOURCE, size=SIZE, md5=MD5,
             filename="db.tar.gz", chunk_size=CHUNK):
    archive = root / filename
    parts = root / "download_parts"
    parts.mkdir(parents=True, exist_ok=True)
    # A first interrupted sequential download can be reused as a partial range.
    if archive.exists() and archive.stat().st_size < size:
        old = archive.stat().st_size
        with archive.open("rb") as src:
            for start in range(0, old, chunk_size):
                part = parts / f"{start:012d}.part"
                n = min(chunk_size, old - start)
                if not part.exists():
                    part.write_bytes(src.read(n))
                else:
                    src.seek(n, 1)
        archive.unlink()
    if not archive.exists():
        ranges = [(i, min(i + chunk_size, size) - 1) for i in range(0, size, chunk_size)]
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(fetch_range, parts, a, b, source, size) for a, b in ranges]
            paths = [job.result() for job in jobs]
        # Check persisted parts again before assembling; final MD5 remains the
        # authoritative content check even when every range had a 206 response.
        for (start, end), path in zip(ranges, paths):
            if path.stat().st_size != end - start + 1:
                fetch_range(parts, start, end, source, size)
        with archive.open("wb") as out:
            for path in paths:
                with path.open("rb") as src:
                    shutil.copyfileobj(src, out)
    if archive.stat().st_size != size or digest(archive, "md5") != md5:
        raise ValueError("Archive did not match official size/MD5")
    print("Archive official MD5 verified", flush=True)
    return archive


def extract(archive, root, loci):
    records = []
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            p = PurePosixPath(member.name)
            if len(p.parts) != 4 or p.parts[:2] != ("db", "loci"):
                continue
            if p.parts[2] not in loci or p.name not in MEMBERS:
                continue
            if not member.isfile():
                raise ValueError(f"Selected member is not regular: {member.name}")
            destination = root.joinpath(*p.parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, destination.open("wb") as out:
                shutil.copyfileobj(src, out)
            records.append({"archive_member": member.name, "bytes": member.size,
                            "sha256": digest(destination)})
    present = {r["archive_member"] for r in records}
    for locus in loci:
        for filename in {"ref.bed", "haplotypes.fa.gz"}:
            if f"db/loci/{locus}/{filename}" not in present:
                raise ValueError(f"Missing required archive member for {locus}: {filename}")
    return records


def fetch_eval_source(root):
    directory = root / "upstream_v0.17.3"
    directory.mkdir(exist_ok=True)
    records = []
    for relpath, sha in EVAL_HASHES.items():
        destination = directory / Path(relpath).name
        url = f"https://raw.githubusercontent.com/tprodanov/locityper/{EVAL_COMMIT}/{relpath}"
        if not destination.exists():
            with urllib.request.urlopen(url, timeout=60) as response:
                destination.write_bytes(response.read())
        if digest(destination) != sha:
            raise ValueError(f"Upstream evaluator hash mismatch: {destination}")
        records.append({"path": relpath, "source_url": url, "sha256": sha})
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=6, choices=range(1, 7))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.cases.open() as inp:
        loci = sorted({row["locus"] for row in csv.DictReader(inp, delimiter="\t")})
    archive = download(args.output, args.workers)
    files = extract(archive, args.output, loci)
    upstream = fetch_eval_source(args.output)
    manifest = {"source_url": SOURCE, "record": "https://zenodo.org/records/15879519",
                "license": "CC-BY-4.0", "official_md5": MD5,
                "bytes": SIZE, "sha256": digest(archive), "selected_loci": loci,
                "files": files,
                "archive_readme_version": "v0.17.3",
                "reference_assembly": "GRCh38",
                "reference_evidence": REFERENCE_EVIDENCE,
                "selection_input": {"path": str(args.cases), "sha256": digest(args.cases)},
                "archive_readme_discrepancy": "README advertises haplotypes.paf.gz, but the checksum-verified archive contains no PAF files. Original archived numeric outcomes cannot be independently replicated until corresponding alignments are recovered or regenerated and validated.",
                "missing_advertised_members": [f"db/loci/{locus}/haplotypes.paf.gz" for locus in loci],
                "evaluation_source": {
                    "repository": "https://github.com/tprodanov/locityper",
                    "tag": "v0.17.3", "commit": EVAL_COMMIT,
                    "files": upstream,
                    "note": "The archive README calls benchmark v0.17.3, while table comments retain v0.17.0 input filenames. Match archived numeric values explicitly."
                },
                "note": "Full archive checksum verified before any extracted panel analysis."}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Extracted {len(files)} selected files for {len(loci)} loci", flush=True)


if __name__ == "__main__":
    main()
