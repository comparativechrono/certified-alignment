#!/usr/bin/env python3
"""Create outcome-blind identity splits for the frozen retrospective study."""
import csv
import gzip
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/panel_adequacy"


def digest(kind, value):
    return hashlib.sha256(f"panel-native-v1|{kind}|{value}".encode()).hexdigest()


def assignments(values, kind):
    values = sorted(set(values), key=lambda x: digest(kind, x))
    a, b = len(values) // 2, len(values) // 4
    return {x: "development" if i < a else "calibration" if i < a+b else "test"
            for i, x in enumerate(values)}


def main():
    with gzip.open(DATA / "locityper/loo_illumina.csv.gz", "rt") as fh:
        rows = csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t")
        # Only identifiers and row type are inspected, never outcomes or diagnostics.
        keys = [(row["sample"], row["locus"]) for row in rows if row["query_type"] == "gt"]
    assert len(keys) == len(set(keys)) == 13040
    donors, loci = {x[0] for x in keys}, {x[1] for x in keys}
    assert len(donors) == 40 and len(loci) == 326
    with open(DATA / "locus_groups.tsv") as fh:
        group_rows = list(csv.DictReader(fh, delimiter="\t"))
    groups = {row["locus"]: row["locus_group"] for row in group_rows}
    assert set(groups) == loci and len(set(groups.values())) == 317
    donor_split = assignments(donors, "donor")
    group_split = assignments(groups.values(), "locus")
    with open(DATA / "donor_splits.tsv", "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "split", "assignment_hash"])
        for donor in sorted(donors):
            writer.writerow([donor, donor_split[donor], digest("donor", donor)])
    with open(DATA / "locus_splits.tsv", "w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["locus", "locus_group", "split", "assignment_hash"])
        for locus in sorted(loci):
            group = groups[locus]
            writer.writerow([locus, group, group_split[group], digest("locus", group)])
    for split in ("development", "calibration", "test"):
        nd = sum(x == split for x in donor_split.values())
        ng = sum(x == split for x in group_split.values())
        nl = sum(group_split[x] == split for x in groups.values())
        print(f"{split}: donors={nd}, locus_groups={ng}, loci={nl}, diagonal_cases={nd*nl}")


if __name__ == "__main__":
    main()
