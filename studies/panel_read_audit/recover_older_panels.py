#!/usr/bin/env python3
"""Recover official older PAFs only if panel sequences and identities match.

The v0.17 archive's missing PAFs motivate this source recovery. A checksum-valid
old archive is not enough: normalized sequences, aliases and bounds are checked
against every selected v0.17 locus before any older alignments can be used.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile

from fetch_panels import download, digest
from resolve_panel_haplotypes import load_locus

SOURCE = "https://zenodo.org/records/12607584/files/db_327.tar.gz?download=1"
SIZE = 934483715
MD5 = "7d5982b60529c8457fd1cc72fcc8ffac"
MEMBERS = {"ref.bed", "haplotypes.fa.gz", "discarded_haplotypes.txt",
           "haplotypes.paf.gz", "haplotypes.paf", "distances.bin"}


def normalized(directory):
    sequences, aliases = load_locus(directory)
    return {"sequence_sha256": {name: hashlib.sha256(seq.upper().encode()).hexdigest()
                                 for name, seq in sorted(sequences.items())},
            "aliases": dict(sorted(aliases.items())),
            "ref_bed": (directory / "ref.bed").read_text().split()}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--current-panel-root", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with args.cases.open() as inp:
        loci = sorted({row["locus"] for row in csv.DictReader(inp, delimiter="\t")})
    archive = download(args.output, 6, source=SOURCE, size=SIZE, md5=MD5,
                       filename="db_327.tar.gz", chunk_size=16 * 1024 * 1024)
    records, paf_members = [], []
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar:
            p = PurePosixPath(member.name)
            if p.name.endswith(".paf") or p.name.endswith(".paf.gz"):
                paf_members.append(member.name)
            if len(p.parts) != 5 or p.parts[:3] != ("dbs", "327", "loci"):
                continue
            locus = p.parts[3]
            if locus not in loci or p.name not in MEMBERS:
                continue
            if not member.isfile():
                raise ValueError(f"Selected member is not regular: {member.name}")
            destination = args.output / "db/loci" / locus / p.name
            destination.parent.mkdir(parents=True, exist_ok=True)
            with tar.extractfile(member) as src, destination.open("wb") as out:
                shutil.copyfileobj(src, out)
            records.append({"archive_member": member.name, "local_member": f"db/loci/{locus}/{p.name}",
                            "bytes": member.size, "sha256": digest(destination)})
    compatibility = []
    for locus in loci:
        old = normalized(args.output / "db/loci" / locus)
        current = normalized(args.current_panel_root / "db/loci" / locus)
        compatibility.append({"locus": locus,
                              "normalized_sequence_ids_and_hashes_match": old["sequence_sha256"] == current["sequence_sha256"],
                              "representatives_and_aliases_match": old["aliases"] == current["aliases"],
                              "reference_bed_matches": old["ref_bed"] == current["ref_bed"],
                              "old": old, "current": current})
    all_match = all(c[k] for c in compatibility for k in (
        "normalized_sequence_ids_and_hashes_match", "representatives_and_aliases_match", "reference_bed_matches"))
    report = {"source_url": SOURCE, "record": "https://zenodo.org/records/12607584",
              "license": "CC-BY-4.0", "official_md5": MD5, "bytes": SIZE,
              "sha256": digest(archive), "selected_loci": loci, "files": records,
              "archive_paf_count": len(paf_members), "all_selected_panels_compatible": all_match,
              "compatibility": compatibility,
              "note": "Panel compatibility is necessary but insufficient. Recovered PAF distances must also reproduce archived full/LOO/available metrics before being treated as original benchmark alignments."}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"archive_paf_count": len(paf_members), "all_selected_panels_compatible": all_match,
                      "selected_file_count": len(records)}), flush=True)


if __name__ == "__main__":
    main()
