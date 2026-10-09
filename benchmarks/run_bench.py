#!/usr/bin/env python3
"""Deterministic synthetic benchmark for the certalign prototype.

This compares modes of one implementation, not external alignment libraries.
Reported per-case times are medians across repetitions. Long-case agreement
does not replace an independently checked optimum.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


MODES = ("full", "linear", "cost", "anchor-linear", "anchor-cost")
SCENARIOS = (
    "uniform_errors", "clustered_errors", "long_insertion", "long_deletion",
    "tandem_repeats", "compensating_indels", "high_similarity",
)
METRICS = (
    "distance", "lower_bound", "rounds", "dp_cells", "certificate_edges",
    "probe_cells", "predictor_ms", "dp_ms", "certificate_ms", "total_ms",
    "peak_cells",
)
REQUIRED = set(METRICS) | {"id", "certified_exact", "cigar", "trace"}
BASE_COLUMNS = (
    "case_id", "mode", "repeat", "batch_index", "scenario", "source_length",
    "seed", "x_length", "y_length", "input_sha256",
)
CIGAR_RE = re.compile(r"([1-9][0-9]*)([MID=X])")


def csv_values(value: str, cast=str):
    try:
        result = [cast(item.strip()) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    if not result:
        raise argparse.ArgumentTypeError("provide at least one comma-separated value")
    if len(result) != len(set(result)):
        raise argparse.ArgumentTypeError("values must be distinct")
    return result


def dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(x: str, rng: random.Random, substitution: float, insertion: float,
           deletion: float, interval=None):
    output = []
    edits = {"substitutions": 0, "inserted_bases": 0, "deleted_bases": 0}
    start, stop = interval if interval is not None else (0, len(x))
    for i, base in enumerate(x):
        if not start <= i < stop:
            output.append(base)
            continue
        if rng.random() < deletion:
            edits["deleted_bases"] += 1
        elif rng.random() < substitution:
            output.append(rng.choice("ACGT".replace(base, "")))
            edits["substitutions"] += 1
        else:
            output.append(base)
        if rng.random() < insertion:
            output.append(rng.choice("ACGT"))
            edits["inserted_bases"] += 1
    return "".join(output), edits


def generate_case(length: int, seed: int, scenario: str):
    label = f"{scenario}-n{length}-s{seed}"
    derived_seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")
    rng = random.Random(derived_seed)
    x = dna(rng, length)
    details = {}
    if scenario == "uniform_errors":
        y, details = mutate(x, rng, .02, .015, .015)
    elif scenario == "clustered_errors":
        start, stop = length * 2 // 5, length * 3 // 5
        y, details = mutate(x, rng, .12, .06, .06, (start, stop))
        details["error_interval"] = [start, stop]
    elif scenario == "long_insertion":
        gap, position = max(16, length // 8), length // 2
        y = x[:position] + dna(rng, gap) + x[position:]
        details = {"inserted_bases": gap, "position": position}
    elif scenario == "long_deletion":
        gap, position = max(16, length // 8), length // 2
        y = x[:position] + x[position + gap:]
        details = {"deleted_bases": gap, "position": position}
    elif scenario == "tandem_repeats":
        motif = dna(rng, rng.randint(2, 6))
        flank = length // 8
        repeat_length = length - 2 * flank
        middle = (motif * math.ceil(repeat_length / len(motif)))[:repeat_length]
        x = x[:flank] + middle + x[-flank:]
        extra = motif * max(2, length // (48 * len(motif)))
        y = x[:flank] + extra + x[flank:]
        y, edits = mutate(y, rng, .004, 0., 0.)
        details = {"motif": motif, "repeat_length": repeat_length,
                   "extra_repeat_bases": len(extra), **edits}
    elif scenario == "compensating_indels":
        gap, deletion_position = max(12, length // 10), length // 4
        shortened = x[:deletion_position] + x[deletion_position + gap:]
        insertion_position = len(shortened) * 3 // 4
        y = shortened[:insertion_position] + dna(rng, gap) + shortened[insertion_position:]
        details = {"deleted_bases": gap, "inserted_bases": gap,
                   "deletion_position": deletion_position,
                   "insertion_position_after_deletion": insertion_position}
    elif scenario == "high_similarity":
        y, details = mutate(x, rng, .002, .001, .001)
    else:
        raise ValueError(f"unknown scenario: {scenario}")
    digest = hashlib.sha256((x + "\n" + y).encode()).hexdigest()
    return {"id": label, "x": x, "y": y, "scenario": scenario,
            "source_length": length, "seed": seed, "derived_seed": derived_seed,
            "x_length": len(x), "y_length": len(y), "input_sha256": digest,
            "generation": details}


def check_cigar(x: str, y: str, cigar: str, expected: int):
    if not isinstance(cigar, str):
        return "CIGAR is not a string"
    matches = list(CIGAR_RE.finditer(cigar))
    if "".join(match.group(0) for match in matches) != cigar:
        return "CIGAR contains invalid or unsupported operations"
    i = j = score = 0
    for match in matches:
        count, op = int(match.group(1)), match.group(2)
        if op in "M=X":
            if i + count > len(x) or j + count > len(y):
                return "CIGAR consumes beyond an input boundary"
            mismatches = sum(a != b for a, b in zip(x[i:i + count], y[j:j + count]))
            if op == "=" and mismatches:
                return "CIGAR equality operation contains a mismatch"
            if op == "X" and mismatches != count:
                return "CIGAR mismatch operation contains an equality"
            score += mismatches
            i, j = i + count, j + count
        elif op == "I":
            j, score = j + count, score + count
        elif op == "D":
            i, score = i + count, score + count
    if (i, j) != (len(x), len(y)):
        return f"CIGAR consumes {(i, j)}, expected {(len(x), len(y))}"
    if score != expected:
        return f"CIGAR score {score} differs from reported distance {expected}"
    return None


def command_for(args, mode):
    return [str(args.binary), "--mode", mode, "--width", str(args.width),
            "--block", str(args.block), "--seed-k", str(args.seed_k),
            "--epsilon", str(args.epsilon)]


def run_batch(args, mode, cases, repeat, batch_index):
    command = command_for(args, mode)
    input_text = "".join(f"{case['id']}\t{case['x']}\t{case['y']}\n" for case in cases)
    started = time.perf_counter()
    try:
        result = subprocess.run(command, input=input_text, text=True, capture_output=True,
                                timeout=args.timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        return [], {"batch_index": batch_index, "mode": mode, "repeat": repeat,
                    "command": command, "case_count": len(cases),
                    "wall_ms": (time.perf_counter() - started) * 1000,
                    "timeout": True}, [f"{mode} repeat {repeat}: timeout after {args.timeout}s"]
    batch = {"batch_index": batch_index, "mode": mode, "repeat": repeat,
             "command": command, "case_count": len(cases),
             "wall_ms": (time.perf_counter() - started) * 1000,
             "returncode": result.returncode, "stderr": result.stderr[-4000:]}
    failures = []
    if result.returncode:
        failures.append(f"{mode} repeat {repeat}: process exited {result.returncode}")
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if len(lines) != len(cases):
        failures.append(f"{mode} repeat {repeat}: expected {len(cases)} JSON rows, got {len(lines)}")
    expected = {case["id"]: case for case in cases}
    seen = set()
    records = []
    for line_number, line in enumerate(lines, 1):
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            failures.append(f"{mode} repeat {repeat} line {line_number}: invalid JSON: {exc}")
            continue
        if not isinstance(payload, dict):
            failures.append(f"{mode} repeat {repeat} line {line_number}: expected JSON object")
            continue
        case_id = payload.get("id")
        if case_id not in expected or case_id in seen:
            failures.append(f"{mode} repeat {repeat}: unexpected or duplicate id {case_id!r}")
            continue
        seen.add(case_id)
        case = expected[case_id]
        prefix = f"{case_id} {mode} repeat {repeat}"
        if "mode" in payload and payload["mode"] != mode:
            failures.append(f"{prefix}: engine reported mode {payload['mode']!r}")
        absent = REQUIRED - payload.keys()
        if absent:
            failures.append(f"{prefix}: missing fields {sorted(absent)}")
        for name in METRICS:
            value = payload.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                failures.append(f"{prefix}: invalid nonnegative metric {name}={value!r}")
        distance = payload.get("distance")
        valid_distance = isinstance(distance, int) and not isinstance(distance, bool) and distance >= 0
        if not valid_distance:
            failures.append(f"{prefix}: distance must be a nonnegative integer")
        if args.epsilon == 0:
            if payload.get("certified_exact") is not True:
                failures.append(f"{prefix}: exact run did not report certified_exact=true")
            if payload.get("lower_bound") != distance:
                failures.append(f"{prefix}: exact run has lower_bound != distance")
        if valid_distance and "cigar" in payload:
            error = check_cigar(case["x"], case["y"], payload["cigar"], distance)
            if error:
                failures.append(f"{prefix}: {error}")
        record = dict(payload)
        record.update({"case_id": case_id, "mode": mode, "repeat": repeat,
                       "batch_index": batch_index, **{key: case[key] for key in
                       ("scenario", "source_length", "seed", "x_length", "y_length", "input_sha256")}})
        records.append(record)
    for case_id in expected.keys() - seen:
        failures.append(f"{mode} repeat {repeat}: missing output for {case_id}")
    return records, batch, failures


def valid_numbers(rows, metric):
    return [row[metric] for row in rows if isinstance(row.get(metric), (int, float))
            and not isinstance(row.get(metric), bool) and math.isfinite(row[metric])]


def summarize(args, cases, records, batches, failures, elapsed_ms):
    by_case_mode = {}
    for row in records:
        by_case_mode.setdefault((row["case_id"], row["mode"]), []).append(row)
    per_case = []
    for (case_id, mode), rows in sorted(by_case_mode.items()):
        item = {"case_id": case_id, "mode": mode, "scenario": rows[0]["scenario"],
                "source_length": rows[0]["source_length"], "seed": rows[0]["seed"],
                "repetitions": len(rows), "metrics": {}}
        for metric in METRICS:
            values = valid_numbers(rows, metric)
            if values:
                item["metrics"][metric] = {"median": statistics.median(values),
                                             "min": min(values), "max": max(values)}
        per_case.append(item)
    for case in cases:
        distances = [(row["mode"], row["repeat"], row["distance"]) for row in records
                     if row["case_id"] == case["id"] and isinstance(row.get("distance"), int)]
        exact_distances = distances if args.epsilon == 0 else [x for x in distances if x[0] == "full"]
        if len({item[2] for item in exact_distances}) > 1:
            failures.append(f"distance mismatch for {case['id']}: {exact_distances}")
        full = [d for mode, _, d in distances if mode == "full"]
        if full and args.epsilon > 0:
            truth = full[0]
            for row in records:
                if row["case_id"] != case["id"]:
                    continue
                distance, lower = row.get("distance"), row.get("lower_bound")
                if isinstance(distance, (int, float)) and isinstance(lower, (int, float)):
                    if not lower <= truth <= distance:
                        failures.append(f"{case['id']} {row['mode']}: bounds exclude full distance {truth}")
    groups = {}
    for item in per_case:
        key = (item["scenario"], item["source_length"], item["mode"])
        groups.setdefault(key, []).append(item)
    aggregates = []
    for (scenario, length, mode), items in sorted(groups.items()):
        aggregate = {"scenario": scenario, "source_length": length, "mode": mode,
                     "case_count": len(items), "metrics": {}}
        for metric in METRICS:
            values = [item["metrics"][metric]["median"] for item in items if metric in item["metrics"]]
            if values:
                aggregate["metrics"][metric] = {
                    "median_of_case_medians": statistics.median(values),
                    "min_case_median": min(values), "max_case_median": max(values),
                    "sum_of_case_medians": sum(values),
                }
        aggregates.append(aggregate)
    lookup = {(item["case_id"], item["mode"]): item for item in per_case}
    comparisons = []
    for baseline, candidate in (("linear", "cost"), ("anchor-linear", "anchor-cost"),
                                ("linear", "anchor-linear"), ("full", "linear"),
                                ("full", "cost"), ("full", "anchor-linear"), ("full", "anchor-cost")):
        common = [case["id"] for case in cases if (case["id"], baseline) in lookup
                  and (case["id"], candidate) in lookup]
        if not common:
            continue
        base = sum(lookup[(key, baseline)]["metrics"].get("total_ms", {}).get("median", 0) for key in common)
        test = sum(lookup[(key, candidate)]["metrics"].get("total_ms", {}).get("median", 0) for key in common)
        comparisons.append({"baseline": baseline, "candidate": candidate,
                            "matched_case_count": len(common), "matched_case_ids": common,
                            "baseline_sum_median_total_ms": base,
                            "candidate_sum_median_total_ms": test,
                            "speedup_ratio": base / test if test > 0 else None})
    full_ids = sorted({row["case_id"] for row in records if row["mode"] == "full"})
    return {"schema_version": 1, "status": "pass" if not failures else "fail",
            "configuration": {key: str(value) if isinstance(value, Path) else value
                              for key, value in vars(args).items()},
            "case_count": len(cases), "result_count": len(records),
            "benchmark_wall_ms": elapsed_ms,
            "validation": {"failures": failures, "full_reference_case_ids": full_ids,
                           "without_full_reference_case_ids": sorted(set(c["id"] for c in cases) - set(full_ids)),
                           "long_case_evidence": "CIGAR validity and internal certificate plus cross-mode agreement; no independent full optimum"},
            "timing_policy": {"primary": "median internal total_ms across repeated identical inputs",
                              "clock": "engine elapsed wall-clock milliseconds, not CPU cycles",
                              "aggregation": "median of per-case medians; comparisons use sums of medians on identical shared cases",
                              "predictor_overhead": "included in engine total_ms; predictor_ms retained separately",
                              "process_overhead": "batch wall_ms includes process startup and I/O and is reported separately",
                              "external_baselines": "none; no Edlib or WFA2 comparison"},
            "batches": batches, "per_case": per_case, "aggregates": aggregates,
            "paired_mode_comparisons": comparisons}


def save_results(output_dir, records, summary):
    (output_dir / "raw_results.json").write_text(json.dumps(records, indent=2) + "\n")
    keys = set().union(*(record.keys() for record in records)) if records else set()
    columns = [name for name in BASE_COLUMNS if name in keys] + sorted(keys - set(BASE_COLUMNS))
    with (output_dir / "raw_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns or list(BASE_COLUMNS))
        writer.writeheader()
        for record in records:
            writer.writerow({key: json.dumps(value, separators=(",", ":")) if isinstance(value, (dict, list)) else value
                             for key, value in record.items()})
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "build" / "certalign")
    parser.add_argument("--lengths", type=lambda s: csv_values(s, int), default=[256, 1024, 4096], help="comma-separated source lengths")
    parser.add_argument("--seeds", type=lambda s: csv_values(s, int), default=[17, 29, 43], help="comma-separated deterministic input seeds")
    parser.add_argument("--repeats", type=int, default=3, help="repeated timing measurements per identical case/mode")
    parser.add_argument("--modes", type=csv_values, default=list(MODES), help="comma-separated modes")
    parser.add_argument("--output-dir", type=Path, default=root / "benchmarks" / "results")
    parser.add_argument("--width", type=int, default=8)
    parser.add_argument("--block", type=int, default=64)
    parser.add_argument("--seed-k", type=int, default=11)
    parser.add_argument("--epsilon", type=float, default=0.)
    parser.add_argument("--full-max-length", type=int, default=1024, help="full DP only if BOTH input lengths are at most this value")
    parser.add_argument("--timeout", type=float, default=90., help="timeout in seconds for one mode/repetition batch")
    parser.add_argument("--generate-only", action="store_true", help="write deterministic inputs and metadata without invoking the engine")
    args = parser.parse_args(argv)
    if any(length < 32 for length in args.lengths):
        parser.error("lengths must be at least 32")
    if (args.repeats < 1 or args.width < 1 or args.block < 1 or args.seed_k < 1
            or not math.isfinite(args.timeout) or args.timeout <= 0
            or not math.isfinite(args.epsilon) or args.epsilon < 0):
        parser.error("repeats, width, block, seed-k and timeout must be positive; epsilon must be nonnegative")
    if not set(args.modes) <= set(MODES):
        parser.error(f"modes must be selected from {','.join(MODES)}")
    args.binary = args.binary.resolve()
    args.output_dir = args.output_dir.resolve()
    if not args.generate_only and not args.binary.is_file():
        parser.error(f"binary does not exist: {args.binary}; build it first or use --generate-only")
    cases = [generate_case(length, seed, scenario) for length in args.lengths
             for seed in args.seeds for scenario in SCENARIOS]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "pairs.tsv").write_text("".join(f"{case['id']}\t{case['x']}\t{case['y']}\n" for case in cases))
    metadata = {"schema_version": 1, "generator": "Python random.Random seeded by SHA256(case id), ACGT alphabet",
                "generation_edits_are_optimal_distances": False,
                "cases": [{key: value for key, value in case.items() if key not in {"x", "y"}} for case in cases]}
    (args.output_dir / "cases.json").write_text(json.dumps(metadata, indent=2) + "\n")
    if args.generate_only:
        print(f"Generated {len(cases)} cases in {args.output_dir}")
        return 0
    records, batches, failures = [], [], []
    binary_sha256 = hashlib.sha256(args.binary.read_bytes()).hexdigest()
    started = time.perf_counter()
    for repeat in range(args.repeats):
        # Rotate mode order to avoid always measuring one method first.
        offset = repeat % len(args.modes)
        mode_order = args.modes[offset:] + args.modes[:offset]
        for mode in mode_order:
            selected = [case for case in cases if mode != "full" or max(case["x_length"], case["y_length"]) <= args.full_max_length]
            if not selected:
                continue
            print(f"repeat {repeat + 1}/{args.repeats}: {mode}, {len(selected)} cases", file=sys.stderr, flush=True)
            try:
                rows, batch, errors = run_batch(args, mode, selected, repeat, len(batches))
            except OSError as exc:
                failures.append(f"failed to execute {args.binary}: {exc}")
                rows, batch, errors = [], {"mode": mode, "repeat": repeat, "execution_error": str(exc)}, []
            records.extend(rows)
            batches.append(batch)
            failures.extend(errors)
            if failures:
                break
        if failures:
            break
    elapsed_ms = (time.perf_counter() - started) * 1000
    if not records and not failures:
        failures.append("no eligible cases were run; check modes and full-max-length")
    summary = summarize(args, cases, records, batches, failures, elapsed_ms)
    summary["engine_binary_sha256"] = binary_sha256
    summary["pairs_tsv_sha256"] = hashlib.sha256((args.output_dir / "pairs.tsv").read_bytes()).hexdigest()
    save_results(args.output_dir, records, summary)
    print(json.dumps({"status": summary["status"], "case_count": len(cases),
                      "result_count": len(records), "wall_seconds": elapsed_ms / 1000,
                      "failure_count": len(failures), "output_dir": str(args.output_dir)}, indent=2))
    for failure in failures[:10]:
        print(f"FAIL: {failure}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
