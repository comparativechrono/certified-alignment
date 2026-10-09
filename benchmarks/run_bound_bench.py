#!/usr/bin/env python3
"""Fixed v2 lower-bound ablation, preserving the v1 inputs and predictors.

The held-out indel cases and controls are fixed before measuring any v2 result.
This compares variants of one implementation, not external alignment libraries.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

from run_bench import METRICS, check_cigar, dna


MODES = ("linear", "cost", "anchor-linear", "anchor-cost")
BOUNDS = ("length", "composition", "qgram")
SETTINGS = tuple((mode, bound) for mode in MODES for bound in BOUNDS)
SEEDS = (101, 103, 107, 109, 113)
INDEL_LENGTHS = (16, 64, 256)
SOURCE_LENGTH = 2048
PREDICTOR_SHA256 = "930a1a18658e3a0b8b651291199afb93912054672daafd7163eb20a112d1156a"
ORIGINAL_PAIRS_SHA256 = "2ef3d5b6c677b87a81500d778a5e0352c2e5c9d1ac3d52c5cd3f43b620ede585"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_original(directory):
    pair_path = directory / "pairs.tsv"
    if sha256(pair_path) != ORIGINAL_PAIRS_SHA256:
        raise ValueError("original pairs.tsv differs from the published v1 inputs")
    metadata = json.loads((directory / "cases.json").read_text())["cases"]
    by_id = {item["id"]: item for item in metadata}
    cases = []
    for line in pair_path.read_text().splitlines():
        case_id, x, y = line.split("\t")
        case = dict(by_id[case_id], x=x, y=y, suite="original")
        if hashlib.sha256((x + "\n" + y).encode()).hexdigest() != case["input_sha256"]:
            raise ValueError(f"original metadata digest mismatch for {case_id}")
        if (len(x), len(y)) != (case["x_length"], case["y_length"]):
            raise ValueError(f"original metadata length mismatch for {case_id}")
        cases.append(case)
    if len(cases) != 63 or len(by_id) != 63:
        raise ValueError("the original suite must contain exactly 63 cases")
    return cases


def make_heldout(seed, scenario, gap):
    label = f"heldout-{scenario}-n{SOURCE_LENGTH}-g{gap}-s{seed}"
    derived_seed = int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")
    rng = random.Random(derived_seed)
    x = dna(rng, SOURCE_LENGTH)
    details = {"indel_length": gap, "substitutions": 0}
    if scenario.startswith("repeat_"):
        # The 1,024-base repeated tract covers the event in both controls.
        motif = dna(rng, 4)
        start, stop = SOURCE_LENGTH // 4, SOURCE_LENGTH * 3 // 4
        x = x[:start] + motif * ((stop - start) // len(motif)) + x[stop:]
        position = SOURCE_LENGTH // 2
        inserted = motif * (gap // len(motif))
        details.update({"motif": motif, "repeat_interval": [start, stop]})
    else:
        # Positions vary reproducibly; none depend on measured algorithm output.
        position = rng.randrange(SOURCE_LENGTH // 3, SOURCE_LENGTH // 2)
        inserted = dna(rng, gap)
    if scenario in {"insertion", "repeat_insertion"}:
        y = x[:position] + inserted + x[position:]
        details.update({"inserted_bases": gap, "deleted_bases": 0, "position": position})
    elif scenario in {"deletion", "repeat_deletion"}:
        y = x[:position] + x[position + gap:]
        details.update({"inserted_bases": 0, "deleted_bases": gap, "position": position})
    elif scenario == "compensating":
        deletion_position = SOURCE_LENGTH // 4 + rng.randrange(64)
        short = x[:deletion_position] + x[deletion_position + gap:]
        insertion_position = len(short) * 3 // 4 + rng.randrange(64)
        y = short[:insertion_position] + inserted + short[insertion_position:]
        details.update({"inserted_bases": gap, "deleted_bases": gap,
                        "deletion_position": deletion_position,
                        "insertion_position_after_deletion": insertion_position})
    else:
        raise ValueError(scenario)
    case = {"id": label, "x": x, "y": y, "suite": "heldout",
            "scenario": scenario, "source_length": SOURCE_LENGTH, "seed": seed,
            "derived_seed": derived_seed, "x_length": len(x), "y_length": len(y),
            "input_sha256": hashlib.sha256((x + "\n" + y).encode()).hexdigest(),
            "generation": details}
    if scenario != "compensating":
        case["analytic_optimal_distance"] = gap
        case["analytic_optimum_reason"] = "length-difference lower bound meets the constructed insertion/deletion path"
    return case


def write_inputs(output_dir, cases):
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "pairs.tsv").write_text("".join(
        f"{case['id']}\t{case['x']}\t{case['y']}\n" for case in cases))
    metadata = {
        "schema_version": 2,
        "original_input_policy": "load the published v1 pairs unchanged; verify SHA256",
        "heldout_generator": "Python random.Random seeded by SHA256(case id), ACGT alphabet",
        "heldout_policy": "fixed seeds, lengths and scenarios before measuring v2; no parameter tuning",
        "generation_edits_are_optimal_distances": False,
        "cases": [{key: value for key, value in case.items() if key not in {"x", "y"}} for case in cases],
    }
    (output_dir / "cases.json").write_text(json.dumps(metadata, indent=2) + "\n")


def command_for(args, mode, bound):
    return [str(args.binary), "--mode", mode, "--bound", bound,
            "--width", str(args.width), "--block", str(args.block),
            "--seed-k", str(args.seed_k), "--epsilon", "0"]


def run_batch(args, setting, cases, repeat, batch_index, phase):
    mode, bound = setting
    command = command_for(args, mode, bound)
    input_text = "".join(f"{c['id']}\t{c['x']}\t{c['y']}\n" for c in cases)
    started = time.perf_counter()
    batch = {"batch_index": batch_index, "mode": mode, "bound": bound,
             "repeat": repeat, "phase": phase, "command": command, "case_count": len(cases)}
    try:
        process = subprocess.run(command, input=input_text, text=True, capture_output=True,
                                 timeout=args.timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        batch.update(wall_ms=(time.perf_counter() - started) * 1000, error=str(exc))
        return [], batch, [f"{mode}/{bound} {phase} {repeat}: {exc}"]
    batch.update(wall_ms=(time.perf_counter() - started) * 1000,
                 returncode=process.returncode, stderr=process.stderr[-4000:])
    failures, records = [], []
    prefix = f"{mode}/{bound} {phase} {repeat}"
    if process.returncode:
        failures.append(f"{prefix}: process exited {process.returncode}")
    lines = [line for line in process.stdout.splitlines() if line.strip()]
    if len(lines) != len(cases):
        failures.append(f"{prefix}: expected {len(cases)} JSON rows; got {len(lines)}")
    expected, seen = {case["id"]: case for case in cases}, set()
    for line_number, line in enumerate(lines, 1):
        try:
            payload = json.loads(line)
        except json.JSONDecodeError as exc:
            failures.append(f"{prefix} line {line_number}: invalid JSON: {exc}")
            continue
        if not isinstance(payload, dict):
            failures.append(f"{prefix} line {line_number}: expected JSON object")
            continue
        case_id = payload.get("id")
        if case_id not in expected or case_id in seen:
            failures.append(f"{prefix}: unexpected or duplicate id {case_id!r}")
            continue
        seen.add(case_id)
        case = expected[case_id]
        label = f"{case_id} {prefix}"
        if payload.get("mode") != mode or payload.get("bound") != bound:
            failures.append(f"{label}: engine reported mismatched mode or bound")
        if payload.get("error"):
            failures.append(f"{label}: engine error {payload['error']}")
        required = set(METRICS) | {"certified_exact", "cigar", "trace"}
        missing = required - payload.keys()
        if missing:
            failures.append(f"{label}: missing fields {sorted(missing)}")
        for name in METRICS:
            value = payload.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                failures.append(f"{label}: invalid nonnegative metric {name}={value!r}")
        distance = payload.get("distance")
        valid_distance = isinstance(distance, int) and not isinstance(distance, bool) and distance >= 0
        if not valid_distance:
            failures.append(f"{label}: distance must be a nonnegative integer")
        if payload.get("certified_exact") is not True or payload.get("lower_bound") != distance:
            failures.append(f"{label}: did not return equal exact bounds")
        if valid_distance:
            error = check_cigar(case["x"], case["y"], payload.get("cigar"), distance)
            if error:
                failures.append(f"{label}: {error}")
            known = case.get("analytic_optimal_distance")
            if known is not None and distance != known:
                failures.append(f"{label}: reported {distance}, analytic optimum is {known}")
        trace = payload.get("trace")
        if not isinstance(trace, list) or not trace:
            failures.append(f"{label}: missing nonempty round trace")
            trace = []
        if valid_distance:
            for index, item in enumerate(trace):
                lower, upper = item.get("lower"), item.get("upper")
                if isinstance(lower, int) and lower >= 0 and lower > distance:
                    failures.append(f"{label}: trace {index} lower exceeds final distance")
                if isinstance(upper, int) and upper >= 0 and upper < distance:
                    failures.append(f"{label}: trace {index} upper below final distance")
        record = dict(payload)
        record.update({"case_id": case_id, "mode": mode, "bound": bound, "phase": phase,
                       "repeat": repeat, "batch_index": batch_index,
                       **{key: case[key] for key in ("suite", "scenario", "source_length", "seed",
                                                    "x_length", "y_length", "input_sha256")}})
        first_optimal = next((i for i, item in enumerate(trace) if item.get("upper") == distance), None)
        if first_optimal is not None:
            post = sum(item["cells"] for item in trace[first_optimal + 1:])
            record["first_optimal_round"] = first_optimal + 1
            record["post_first_optimal_dp_cells"] = post
            record["post_first_optimal_dp_fraction"] = post / payload["dp_cells"] if payload.get("dp_cells") else 0.
        records.append(record)
    for case_id in expected.keys() - seen:
        failures.append(f"{prefix}: missing output for {case_id}")
    return records, batch, failures


def numeric_keys(records):
    excluded = {"repeat", "batch_index", "source_length", "seed", "n", "m", "x_length", "y_length"}
    return sorted({key for row in records for key, value in row.items()
                   if key not in excluded and isinstance(value, (int, float))
                   and not isinstance(value, bool) and math.isfinite(value)})


def summarize(args, cases, records, batches, failures, elapsed_ms):
    timed = [row for row in records if row["phase"] == "timed"]
    metrics = numeric_keys(timed)
    groups = {}
    for row in timed:
        groups.setdefault((row["case_id"], row["mode"], row["bound"]), []).append(row)
    lookup, per_case = {}, []
    for (case_id, mode, bound), rows in sorted(groups.items()):
        item = {"case_id": case_id, "mode": mode, "bound": bound,
                "suite": rows[0]["suite"], "scenario": rows[0]["scenario"],
                "source_length": rows[0]["source_length"], "seed": rows[0]["seed"],
                "repetitions": len(rows), "metrics": {}}
        for metric in metrics:
            values = [row[metric] for row in rows if isinstance(row.get(metric), (int, float))
                      and not isinstance(row.get(metric), bool) and math.isfinite(row[metric])]
            if values:
                item["metrics"][metric] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
        if len(rows) != args.repeats:
            failures.append(f"{case_id} {mode}/{bound}: expected {args.repeats} repeats; got {len(rows)}")
        lookup[(case_id, mode, bound)] = item
        per_case.append(item)
    expected_count = len(cases) * len(SETTINGS) * args.repeats
    if len(timed) != expected_count:
        failures.append(f"expected {expected_count} timed outputs; got {len(timed)}")
    for case in cases:
        rows = [row for row in records if row["case_id"] == case["id"]]
        if len({row.get("distance") for row in rows}) != 1:
            failures.append(f"cross-setting/repetition distance mismatch for {case['id']}")

    # Every scope uses identical cases for paired comparisons and sums case medians.
    scopes = [("all", cases)]
    for suite in ("original", "heldout"):
        selected = [case for case in cases if case["suite"] == suite]
        scopes.append((suite, selected))
        for scenario in sorted({case["scenario"] for case in selected}):
            scenario_cases = [case for case in selected if case["scenario"] == scenario]
            scopes.append((f"{suite}/{scenario}", scenario_cases))
            if suite == "heldout":
                for gap in sorted({case["generation"]["indel_length"] for case in scenario_cases}):
                    scopes.append((f"{suite}/{scenario}/g{gap}",
                                   [case for case in scenario_cases if case["generation"]["indel_length"] == gap]))
    aggregates, bound_comparisons, predictor_comparisons = [], [], []
    for scope, selected in scopes:
        sums = {}
        for mode, bound in SETTINGS:
            items = [lookup[(case["id"], mode, bound)] for case in selected
                     if (case["id"], mode, bound) in lookup]
            aggregated = {"scope": scope, "mode": mode, "bound": bound,
                          "case_count": len(items), "metrics": {}}
            for metric in metrics:
                values = [item["metrics"][metric]["median"] for item in items if metric in item["metrics"]]
                if values:
                    aggregated["metrics"][metric] = {"sum_of_case_medians": sum(values),
                                                      "median_of_case_medians": statistics.median(values)}
            # This ratio weights by work; the mean of per-case fractions does not.
            total_cells = aggregated["metrics"].get("dp_cells", {}).get("sum_of_case_medians", 0)
            post_cells = aggregated["metrics"].get("post_first_optimal_dp_cells", {}).get("sum_of_case_medians", 0)
            aggregated["post_first_optimal_fraction_of_total_cells"] = post_cells / total_cells if total_cells else None
            sums[(mode, bound)] = aggregated
            aggregates.append(aggregated)

        def compare(baseline, candidate):
            left, right = sums[baseline], sums[candidate]
            result = {"scope": scope, "baseline_mode": baseline[0], "baseline_bound": baseline[1],
                      "candidate_mode": candidate[0], "candidate_bound": candidate[1],
                      "matched_case_count": len(selected), "metrics": {}}
            for metric in metrics:
                a = left["metrics"].get(metric, {}).get("sum_of_case_medians")
                b = right["metrics"].get(metric, {}).get("sum_of_case_medians")
                if a is not None and b is not None:
                    result["metrics"][metric] = {"baseline_sum_median": a, "candidate_sum_median": b,
                                                 "candidate_over_baseline": b / a if a else None,
                                                 "relative_change": (b - a) / a if a else None}
            a = result["metrics"].get("total_ms", {}).get("baseline_sum_median", 0)
            b = result["metrics"].get("total_ms", {}).get("candidate_sum_median", 0)
            result["speedup_ratio"] = a / b if b else None
            return result

        for mode in MODES:
            for bound in BOUNDS[1:]:
                bound_comparisons.append(compare((mode, "length"), (mode, bound)))
        for bound in BOUNDS:
            for baseline, candidate in (("linear", "cost"), ("anchor-linear", "anchor-cost")):
                predictor_comparisons.append(compare((baseline, bound), (candidate, bound)))

    return {"schema_version": 2, "status": "fail" if failures else "pass",
            "configuration": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "fixed_configuration": {"modes": MODES, "bounds": BOUNDS, "epsilon": 0,
                                    "heldout_source_length": SOURCE_LENGTH, "heldout_seeds": SEEDS,
                                    "heldout_indel_lengths": INDEL_LENGTHS, "warmups_per_setting": 1},
            "case_count": len(cases), "timed_result_count": len(timed), "result_count_including_warmup": len(records),
            "benchmark_wall_ms": elapsed_ms, "validation": {
                "failures": failures, "failure_count": len(failures),
                "evidence": "All CIGARs replay; exact certificates equal final scores; cross-setting/repetition agreement. Forty held-out cases additionally have an analytic optimum from a pure insertion/deletion construction. No independent full-DP oracle is run in this harness.",
                "analytic_reference_case_count": sum("analytic_optimal_distance" in case for case in cases),
                "trace_check": "Each finite trace lower/upper brackets the final certified distance; this is not independent proof of optimum.",
            },
            "timing_policy": {
                "primary": "per-case median internal total_ms over five repetitions by default",
                "aggregation": "paired sums of per-case medians on identical cases, separately by suite/scenario",
                "order": "12 mode/bound settings; rotate by 5*repeat positions, reverse odd repetitions",
                "warmup": "one complete untimed batch per setting; retained in raw outputs, excluded from medians",
                "clock": "engine elapsed wall-clock milliseconds, not CPU cycles",
                "process": "fresh process per batch; batch wall_ms separately includes process startup and I/O",
                "overheads": "all engine preprocessing, prediction and certification included in total_ms",
                "external_baselines": "none; no external-library or state-of-art speed claim",
                "post_first_optimal_dp_cells": "cells in later rounds after the first round whose upper equals the final exact score, identified retrospectively; excludes that first successful round",
            },
            "metric_names": metrics, "batches": batches, "per_case": per_case,
            "aggregates": aggregates, "paired_bound_comparisons": bound_comparisons,
            "paired_predictor_comparisons": predictor_comparisons}


def save_results(output_dir, records, summary):
    (output_dir / "raw_results.json").write_text(json.dumps(records, indent=2) + "\n")
    keys = set().union(*(row.keys() for row in records)) if records else set()
    leading = [name for name in ("case_id", "suite", "scenario", "mode", "bound", "phase", "repeat", "batch_index") if name in keys]
    columns = leading + sorted(keys - set(leading))
    with (output_dir / "raw_results.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in records:
            writer.writerow({key: json.dumps(value, separators=(",", ":")) if isinstance(value, (dict, list)) else value
                             for key, value in row.items()})
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "build" / "certalign")
    parser.add_argument("--original-dir", type=Path, default=root / "results" / "benchmark")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "bound_benchmark")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--width", type=int, default=8)
    parser.add_argument("--block", type=int, default=64)
    parser.add_argument("--seed-k", type=int, default=11)
    parser.add_argument("--timeout", type=float, default=90.)
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args(argv)
    if min(args.repeats, args.width, args.block, args.seed_k) < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("repeats, width, block, seed-k and timeout must be positive")
    args.binary, args.original_dir, args.output_dir = args.binary.resolve(), args.original_dir.resolve(), args.output_dir.resolve()
    if args.original_dir == args.output_dir:
        parser.error("output-dir must differ from original-dir")
    if not args.generate_only and not args.binary.is_file():
        parser.error(f"binary does not exist: {args.binary}")
    predictor_sha = sha256(root / "include" / "predictor.hpp")
    if predictor_sha != PREDICTOR_SHA256:
        parser.error("predictor.hpp differs from the fixed v1 predictor; this experiment requires it unchanged")
    cases = load_original(args.original_dir)
    cases += [make_heldout(seed, scenario, gap) for seed in SEEDS
              for scenario in ("insertion", "deletion", "compensating") for gap in INDEL_LENGTHS]
    cases += [make_heldout(seed, scenario, 64) for seed in SEEDS
              for scenario in ("repeat_insertion", "repeat_deletion")]
    if len(cases) != 118 or len({case["id"] for case in cases}) != len(cases):
        raise ValueError("invalid fixed suite")
    write_inputs(args.output_dir, cases)
    if args.generate_only:
        print(json.dumps({"case_count": len(cases), "original": 63, "heldout": 55,
                          "pairs_sha256": sha256(args.output_dir / "pairs.tsv"),
                          "output_dir": str(args.output_dir)}, indent=2))
        return 0
    binary_sha = sha256(args.binary)
    records, batches, failures = [], [], []
    started = time.perf_counter()
    orders = [("warmup", -1, list(SETTINGS))]
    for repeat in range(args.repeats):
        offset = 5 * repeat % len(SETTINGS)
        order = list(SETTINGS[offset:] + SETTINGS[:offset])
        if repeat % 2:
            order.reverse()
        orders.append(("timed", repeat, order))
    for phase, repeat, order in orders:
        for setting in order:
            print(f"{phase} {repeat + 1}: {'/'.join(setting)}, {len(cases)} cases", file=sys.stderr, flush=True)
            rows, batch, errors = run_batch(args, setting, cases, repeat, len(batches), phase)
            records.extend(rows)
            batches.append(batch)
            failures.extend(errors)
            if failures:
                break
        if failures:
            break
    elapsed_ms = (time.perf_counter() - started) * 1000
    summary = summarize(args, cases, records, batches, failures, elapsed_ms)
    summary.update(engine_binary_sha256=binary_sha, predictor_sha256=predictor_sha,
                   original_pairs_sha256=ORIGINAL_PAIRS_SHA256,
                   pairs_tsv_sha256=sha256(args.output_dir / "pairs.tsv"))
    save_results(args.output_dir, records, summary)
    print(json.dumps({"status": summary["status"], "case_count": len(cases),
                      "timed_result_count": summary["timed_result_count"], "wall_seconds": elapsed_ms / 1000,
                      "failure_count": len(failures), "output_dir": str(args.output_dir)}, indent=2))
    for failure in failures[:10]:
        print(f"FAIL: {failure}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
