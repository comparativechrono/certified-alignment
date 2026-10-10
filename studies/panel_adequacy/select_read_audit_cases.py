#!/usr/bin/env python3
"""Reproduce truth-selected exploratory read-audit cases and deterministic controls.

This is post-hoc case selection from known baseline errors, not a new held-out
performance test. Source truth is used for selection only. No read outcomes are
examined by this script. Coordinates are a grouping BED aid; downstream read
extraction must verify the locus database reference and exact bounds.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/panel_adequacy_baseline"
DATA = ROOT / "data/panel_adequacy"
OUT = ROOT / "data/panel_read_audit"
THRESHOLD = 0.005
ACCEPT = "accept_matched_0.95:unexpl_fraction"
CONTROL_SALT = "panel-read-audit-v1|control"


def read(path: Path, sep: str = "\t") -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader((r for r in handle if not r.startswith("#")), delimiter=sep))


def key(row):
    return row["sample"], row["locus"]


def unique_index(rows):
    result = {key(row): row for row in rows}
    assert len(result) == len(rows), "Duplicate sample/locus key"
    return result


def number(row, column):
    return float(row[column])


def scorable(row):
    return (math.isfinite(number(row, "div")) and number(row, "div") >= 0
        and math.isfinite(number(row, "size")) and number(row, "size") > 0)


def below(row):
    return scorable(row) and number(row, "div") <= THRESHOLD


def write(path, rows):
    assert rows
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def control_hash(case_key, control_key):
    payload = f"{CONTROL_SALT}|{case_key[0]}|{case_key[1]}|{control_key[0]}"
    return hashlib.sha256(payload.encode()).hexdigest()


def main():
    inputs = [
        BASE / "heldout_predictions.csv.gz", BASE / "heldout_truth.csv.gz",
        DATA / "locityper/loo_illumina.csv.gz", DATA / "locityper/full_illumina.csv.gz",
        DATA / "locus_groups.tsv", DATA / "donor_splits.tsv", DATA / "locus_splits.tsv",
    ]
    pred = unique_index([r for r in read(inputs[0], ",") if r["condition"] == "loo"])
    truth = unique_index([r for r in read(inputs[1], ",") if r["condition"] == "loo"])
    loo = unique_index([r for r in read(inputs[2]) if r["query_type"] == "gt"])
    full = unique_index([r for r in read(inputs[3]) if r["query_type"] == "gt"])
    groups = {r["locus"]: r for r in read(inputs[4])}
    donors = {r["sample"]: r["split"] for r in read(inputs[5])}
    loci = {r["locus"]: r["split"] for r in read(inputs[6])}
    assert set(pred) == set(truth) and len(pred) == 806
    assert all(donors[s] == loci[l] == "test" for s, l in pred)
    for k in pred:
        assert number(truth[k], "div") == number(loo[k], "div")
        assert number(truth[k], "avail_div") == number(loo[k], "avail_div")
        for column in ["quality", "total_reads", "unexpl_reads", "weight_dist"]:
            assert number(pred[k], column) == number(loo[k], column)
    selected = sorted((k for k in pred
        if int(pred[k][ACCEPT]) == 1
        and number(truth[k], "div") > THRESHOLD
        and number(truth[k], "avail_div") > THRESHOLD
        and below(full[k])), key=lambda k: (k[1], k[0]))
    assert len(selected) == 22, "Committed baseline no longer produces exactly 22 candidates"
    controls = [k for k in pred if below(loo[k]) and below(full[k])]
    uses = Counter()
    cases, matches = [], []

    def describe(k, role):
        s, l = k
        row = dict(sample=s, locus=l, role=role,
            truth_status="posthoc_sequence_proxy_not_validated_structural_genotype",
            grouping_bed_chrom=groups[l]["chrom"], grouping_bed_start=groups[l]["start"],
            grouping_bed_end=groups[l]["end"], grouping_bed_source="data/panel_adequacy/cosigt/cmrgs_original.bed",
            coordinate_status="grouping_aid_verify_database_reference_before_read_extraction",
            locus_group=groups[l]["locus_group"],
            expected_truth_haplotype_ids=f"{s}.1,{s}.2",
            truth_haplotype_resolution="requires_database_discarded_haplotype_alias_check",
            baseline_fraction_accepted=int(pred[k][ACCEPT]))
        for condition, source in [("loo", loo[k]), ("full", full[k])]:
            for column in ["genotype", "quality", "total_reads", "unexpl_reads", "weight_dist", "warnings", "edit", "size", "div", "qv"]:
                row[f"{condition}_{column}"] = source[column]
            row[f"{condition}_unexpl_fraction"] = number(source, "unexpl_reads") / number(source, "total_reads")
        for column in ["avail_edit", "avail_size", "avail_div", "avail_qv"]:
            row[f"loo_{column}"] = loo[k][column]
        return row

    control_rows = {}
    for k in selected:
        case_id = f"{k[0]}__{k[1]}"
        candidates = [c for c in controls if c[1] == k[1] and c[0] != k[0]]
        match = dict(case_id=case_id, case_sample=k[0], locus=k[1],
            eligible_same_locus_controls=len(candidates), control_sample="", control_id="",
            control_prior_use_count="", control_selection_hash="", match_status="unmatched_no_eligible_heldout_control")
        if candidates:
            c = min(candidates, key=lambda candidate: (uses[candidate], control_hash(k, candidate)))
            match.update(control_sample=c[0], control_id=f"{c[0]}__{c[1]}",
                control_prior_use_count=uses[c], control_selection_hash=control_hash(k, c),
                match_status="matched_reused_control" if uses[c] else "matched_unique_control")
            uses[c] += 1
            control_rows[c] = describe(c, "represented_sequence_proxy_control")
        row = dict(case_id=case_id, **describe(k, "panel_adequacy_candidate"))
        row["control_id"] = match["control_id"]
        row["match_status"] = match["match_status"]
        cases.append(row)
        matches.append(match)
    OUT.mkdir(parents=True, exist_ok=True)
    write(OUT / "cases.tsv", cases)
    write(OUT / "control_manifest.tsv", matches)
    control_details = [dict(control_id=f"{k[0]}__{k[1]}", **control_rows[k], case_match_count=uses[k])
        for k in sorted(control_rows, key=lambda k: (k[1], k[0]))]
    write(OUT / "controls.tsv", control_details)
    manifest = dict(schema_version=1, design="posthoc_exploratory_case_control_audit",
        historical_baseline_commit="9b33cd9",
        selection=f"LOO heldout {ACCEPT}=1 and div>{THRESHOLD} and avail_div>{THRESHOLD} and paired full div<={THRESHOLD}",
        controls=f"Same-locus different-donor heldout case with scorable full and LOO div<={THRESHOLD}; no acceptance requirement; no read-outcome or native-metric ranking",
        matching="Process candidates by (locus,sample); choose minimum prior control reuse count, then SHA256(control salt|case sample|locus|control sample); do not expand beyond heldout cases when unmatched",
        control_hash_salt=CONTROL_SALT,
        cases=len(cases), loci=len({k[1] for k in selected}), donors=len({k[0] for k in selected}),
        matched_cases=sum(bool(m["control_id"]) for m in matches), unique_control_cases=len(control_rows),
        reused_control_matches=sum(m["match_status"] == "matched_reused_control" for m in matches),
        unmatched_case_ids=[m["case_id"] for m in matches if not m["control_id"]],
        warnings=["Truth-selected cases cannot estimate prospective accuracy or demonstrate calibrated uncertainty.",
            "The complete set contains 22 known sequence discrepancies, not 22 proven absent structural alleles.",
            "Reused controls are not independent observations.",
            "Expected sample.1/sample.2 IDs require alias verification against the archived database.",
            "Grouping BED coordinates are not authorized as extraction coordinates without reference verification."],
        inputs=[dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size) for p in inputs],
        outputs=[dict(path=str(p.relative_to(ROOT)), sha256=hashlib.sha256(p.read_bytes()).hexdigest(), bytes=p.stat().st_size)
            for p in [OUT / "cases.tsv", OUT / "control_manifest.tsv", OUT / "controls.tsv"]])
    (OUT / "SELECTION_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ["cases", "loci", "donors", "matched_cases", "unique_control_cases", "reused_control_matches", "unmatched_case_ids"]}, indent=2))


if __name__ == "__main__":
    main()
