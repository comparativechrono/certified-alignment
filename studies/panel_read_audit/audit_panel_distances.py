#!/usr/bin/env python3
"""Independently check archived distances using the unmodified upstream code.

Run only after the read-audit protocol/code freeze. This is an input-compatibility
oracle, not a newly developed measure of correctness or read support.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

from fetch_panels import EVAL_HASHES, EVAL_COMMIT, digest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel-root", type=Path, required=True)
    ap.add_argument("--upstream-root", type=Path,
                    help="Optional pinned evaluator directory when alignments are recovered elsewhere.")
    ap.add_argument("--cases", type=Path, required=True)
    ap.add_argument("--controls", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    upstream = args.upstream_root or args.panel_root / "upstream_v0.17.3"
    for path, expected in EVAL_HASHES.items():
        assert digest(upstream / Path(path).name) == expected, path
    sys.path.insert(0, str(upstream.resolve()))
    from gt_dist import Distances

    rows = []
    for table in (args.cases, args.controls):
        with table.open() as inp:
            rows.extend(csv.DictReader(inp, delimiter="\t"))
    assert len({(r["sample"], r["locus"]) for r in rows}) == len(rows)
    loaded = {}
    results = []
    for row in rows:
        locus = row["locus"]
        directory = args.panel_root / "db/loci" / locus
        if locus not in loaded:
            loaded[locus] = Distances(str(directory / "discarded_haplotypes.txt"),
                                     str(directory / "haplotypes.paf.gz"))
        distances = loaded[locus]
        truth = list(distances.get_sample_haplotypes(row["sample"]))
        assert len(truth) == 2, (row["sample"], locus, truth)
        closest, oracle = distances.find_closest(truth, loo=True)
        full = distances.calc_distance(truth, row["full_genotype"].split(","))
        loo = distances.calc_distance(truth, row["loo_genotype"].split(","))
        checks = {}
        for name, result, prefix in (("full", full, "full_"), ("loo", loo, "loo_"),
                                    ("available", oracle, "loo_avail_")):
            out = list(result.iter_strs())[-1].split("\t")
            actual = dict(zip(("edit", "size", "div", "qv"), out))
            status = {key: float(actual[key]) == float(row[prefix + key])
                      for key in ("edit", "size", "div", "qv")}
            checks[name] = {"computed": actual, "matches_archive": status}
        results.append({"sample": row["sample"], "locus": locus,
                        "truth_haplotype_ids": truth,
                        "closest_available_haplotype_ids": closest,
                        "comparisons": checks})
    mismatches = [{"sample": r["sample"], "locus": r["locus"], "comparison": n,
                   "field": k} for r in results for n, c in r["comparisons"].items()
                  for k, ok in c["matches_archive"].items() if not ok]
    report = {"method": "unmodified upstream Locityper gt_dist.py",
              "upstream_commit": EVAL_COMMIT, "cases": len(results),
              "comparisons": 3 * len(results), "mismatches": mismatches,
              "rows": results}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("cases", "comparisons", "mismatches")}))


if __name__ == "__main__":
    main()
