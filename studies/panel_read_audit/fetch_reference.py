#!/usr/bin/env python3
"""Fetch required GRCh38 chromosomes and verify official sequence MD5 values.

Byte ranges come from the official reference FAI. Complete chromosomes, rather
than padded locus slices, permit ordinary CRAM reference MD5 validation.
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
from pathlib import Path

from fetch_panels import fetch_range, digest

SOURCE = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
SIZE = 3263683042


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--dictionary", type=Path, required=True)
    ap.add_argument("--panel-root", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=6, choices=range(1, 7))
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    parts_dir = args.output_dir / "download_parts"
    parts_dir.mkdir(exist_ok=True)
    with args.cases.open() as inp:
        loci = {row["locus"] for row in csv.DictReader(inp, delimiter="\t")}
    chromosomes = { (args.panel_root / "db/loci" / locus / "ref.bed").read_text().split()[0]
                    for locus in loci }
    dictionary = {}
    for line in args.dictionary.read_text().splitlines():
        if line.startswith("@SQ\t"):
            fields = dict(item.split(":", 1) for item in line.split("\t")[1:])
            dictionary[fields["SN"]] = fields
    reference = []
    for line in args.fai.read_text().splitlines():
        name, n, offset, bases, width = line.split("\t")
        if name in chromosomes:
            n, offset, bases, width = map(int, (n, offset, bases, width))
            assert int(dictionary[name]["LN"]) == n
            end = offset + (n - 1) // bases * width + (n - 1) % bases
            ranges = [(start, min(start + 16 * 1024 * 1024 - 1, end))
                      for start in range(offset, end + 1, 16 * 1024 * 1024)]
            reference.append({"name": name, "length": n, "offset": offset,
                              "line_bases": bases, "line_width": width,
                              "end_inclusive": end, "ranges": ranges,
                              "expected_md5": dictionary[name]["M5"]})
    assert len(reference) == len(chromosomes)
    jobs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for seq in reference:
            for start, end in seq["ranges"]:
                jobs.append(pool.submit(fetch_range, parts_dir, start, end, SOURCE, SIZE))
        for job in concurrent.futures.as_completed(jobs):
            job.result()
    for seq in reference:
        for start, end in seq["ranges"]:
            part = parts_dir / f"{start:012d}.part"
            if part.stat().st_size != end - start + 1:
                fetch_range(parts_dir, start, end, SOURCE, SIZE)
    output = args.output_dir / "selected_GRCh38.fa"
    fai_lines = []
    with output.open("wb") as out:
        for seq in reference:
            out.write((">" + seq["name"] + "\n").encode())
            new_offset = out.tell()
            h, n = hashlib.md5(), 0
            for start, end in seq["ranges"]:
                part = parts_dir / f"{start:012d}.part"
                with part.open("rb") as inp:
                    for block in iter(lambda: inp.read(4 * 1024 * 1024), b""):
                        out.write(block)
                        sequence = block.replace(b"\n", b"").replace(b"\r", b"").upper()
                        h.update(sequence)
                        n += len(sequence)
            out.write(b"\n")
            assert n == seq["length"], (seq["name"], n, seq["length"])
            assert h.hexdigest() == seq["expected_md5"], (seq["name"], h.hexdigest())
            seq["verified_md5"] = h.hexdigest()
            fai_lines.append(f"{seq['name']}\t{n}\t{new_offset}\t{seq['line_bases']}\t{seq['line_width']}")
            print(f"Verified chromosome {seq['name']} MD5", flush=True)
    Path(str(output) + ".fai").write_text("\n".join(fai_lines) + "\n")
    report = {"source_url": SOURCE, "source_bytes": SIZE,
              "fai_source_url": SOURCE + ".fai", "fai_sha256": digest(args.fai),
              "dictionary_source_url": SOURCE.removesuffix(".fa") + ".dict",
              "dictionary_sha256": digest(args.dictionary),
              "selected_chromosomes": reference, "output_bytes": output.stat().st_size,
              "output_sha256": digest(output),
              "output_fai_sha256": digest(Path(str(output) + ".fai")),
              "validation": "All selected complete chromosome sequences match official reference.dict M5 checksums."}
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Validated reference ready: {output}", flush=True)


if __name__ == "__main__":
    main()
