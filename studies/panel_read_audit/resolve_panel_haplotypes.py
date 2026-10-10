#!/usr/bin/env python3
"""Resolve frozen case/control haplotype IDs against checksum-verified FASTAs.

This source compatibility check reports metadata only, without comparing the
haplotype sequences or computing structural/read outcomes.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def load_locus(directory):
    sequences = {}
    name, chunks = None, []
    with gzip.open(directory / "haplotypes.fa.gz", "rt") as inp:
        for line in inp:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        if name is not None:
            sequences[name] = "".join(chunks)
    alias_to_representative = {name: name for name in sequences}
    aliases = directory / "discarded_haplotypes.txt"
    if aliases.exists():
        for line in aliases.read_text().splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            representative, others = line.split("=")
            representative = representative.strip()
            assert representative in sequences
            for name in others.split(","):
                name = name.strip()
                assert name not in alias_to_representative
                alias_to_representative[name] = representative
    return sequences, alias_to_representative


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel-root", type=Path, required=True)
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--controls", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    rows = []
    for table in (args.cases, args.controls):
        with table.open() as inp:
            rows.extend(csv.DictReader(inp, delimiter="\t"))
    loaded, results, missing = {}, [], []
    for row in rows:
        locus, sample = row["locus"], row["sample"]
        if locus not in loaded:
            loaded[locus] = load_locus(args.panel_root / "db/loci" / locus)
        sequences, aliases = loaded[locus]
        groups = {"truth": [sample + ".1", sample + ".2"],
                  "full_prediction": row["full_genotype"].split(","),
                  "loo_prediction": row["loo_genotype"].split(",")}
        resolved = {}
        for label, names in groups.items():
            resolved[label] = []
            for name in names:
                representative = aliases.get(name)
                if representative is None:
                    missing.append({"sample": sample, "locus": locus, "id": name, "role": label})
                    resolved[label].append({"id": name, "status": "missing"})
                    continue
                seq = sequences[representative]
                resolved[label].append({"id": name, "representative": representative,
                                        "alias": name != representative, "length": len(seq),
                                        "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest()})
        results.append({"sample": sample, "locus": locus, "resolved": resolved})
    report = {"rows": results, "case_control_count": len(rows), "missing": missing,
              "note": "Metadata compatibility only. No PAFs occur in the pinned archive, so numeric distance compatibility is not established by this check."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"case_control_count": len(rows), "missing": missing}))


if __name__ == "__main__":
    main()
