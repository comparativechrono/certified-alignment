#!/usr/bin/env python3
"""Fixed v3 ordered-bound ablation on the published 118 pairs plus 30 new pairs.

Reuse v2 CIGAR/certificate validation and raw-result serialization, without
changing its harness or predictor. Inputs and parameters are fixed before
timing; this compares variants of one engine, not external aligners.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import random
import statistics
import sys
import time
from pathlib import Path

import run_bound_bench as previous

MODES = ("linear", "cost", "anchor-linear", "anchor-cost")
BOUNDS = ("length", "qgram", "ordered")
SETTINGS = tuple((mode, bound) for mode in MODES for bound in BOUNDS)
PRIMARY_SHA256 = "f22ff60cf885ebf420da187d86bc33baac3c80dffd0c9e822bd7c6f71b1382b5"
FRESH_SEEDS = (211, 223, 227, 229, 233)
FRESH_GAPS = (32, 128, 384)
FRESH_LENGTH = 2048
ORDERED_METRICS = ("ordered_cells", "ordered_bytes", "ordered_preprocess_ms", "ordered_calls")


def load_primary(directory, original_directory):
    original = previous.load_original(original_directory)
    pairs_path = directory / "pairs.tsv"
    if previous.sha256(pairs_path) != PRIMARY_SHA256:
        raise ValueError("primary 118 pairs differ from the published v2 inputs")
    metadata = json.loads((directory / "cases.json").read_text())["cases"]
    by_id = {case["id"]: case for case in metadata}
    cases = []
    for line in pairs_path.read_text().splitlines():
        case_id, x, y = line.split("\t")
        case = dict(by_id[case_id], x=x, y=y)
        if hashlib.sha256((x + "\n" + y).encode()).hexdigest() != case["input_sha256"]:
            raise ValueError(f"metadata digest mismatch for {case_id}")
        if (len(x), len(y)) != (case["x_length"], case["y_length"]):
            raise ValueError(f"metadata length mismatch for {case_id}")
        cases.append(case)
    if len(cases) != 118 or len(by_id) != 118:
        raise ValueError("primary suite must contain exactly 118 cases")
    existing = {case["id"]: (case["x"], case["y"]) for case in cases}
    if any(existing.get(case["id"]) != (case["x"], case["y"]) for case in original):
        raise ValueError("primary suite does not contain the original 63 unchanged pairs")
    if sum("analytic_optimal_distance" in case for case in cases) != 40:
        raise ValueError("primary suite must preserve its forty analytic references")
    return cases


def make_fresh(seed, gap, repeat_background):
    scenario = "repeat_compensating" if repeat_background else "compensating"
    case_id = f"fresh-v3-{scenario}-n{FRESH_LENGTH}-g{gap}-s{seed}"
    derived_seed = int.from_bytes(hashlib.sha256(case_id.encode()).digest()[:8], "big")
    rng = random.Random(derived_seed)
    x = previous.dna(rng, FRESH_LENGTH)
    details = {"indel_length": gap, "inserted_bases": gap, "deleted_bases": gap,
               "substitutions": 0, "repeat_background": repeat_background}
    if repeat_background:
        motif = previous.dna(rng, 4)
        start, stop = FRESH_LENGTH // 4, FRESH_LENGTH * 3 // 4
        x = x[:start] + motif * ((stop - start) // 4) + x[stop:]
        details.update(motif=motif, repeat_interval=[start, stop])
    deletion_position = FRESH_LENGTH // 4 + rng.randrange(64)
    short = x[:deletion_position] + x[deletion_position + gap:]
    insertion_position = len(short) * 3 // 4 + rng.randrange(64)
    inserted = previous.dna(rng, gap)
    y = short[:insertion_position] + inserted + short[insertion_position:]
    details.update(deletion_position=deletion_position,
                   insertion_position_after_deletion=insertion_position)
    return {"id": case_id, "x": x, "y": y, "suite": "fresh_v3", "scenario": scenario,
            "source_length": FRESH_LENGTH, "seed": seed, "derived_seed": derived_seed,
            "x_length": len(x), "y_length": len(y), "generation": details,
            "input_sha256": hashlib.sha256((x + "\n" + y).encode()).hexdigest()}


def write_inputs(directory, cases):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "pairs.tsv").write_text("".join(
        f"{case['id']}\t{case['x']}\t{case['y']}\n" for case in cases))
    metadata = {"schema_version": 3,
                "primary_policy": "published 118 v2 pairs loaded unchanged and SHA256 verified",
                "original_policy": "original 63 v1 pairs additionally verified by their published SHA256",
                "fresh_policy": "30 new compensated-indel pairs; fixed seeds, gaps, backgrounds and parameters before v3 timing",
                "fresh_generator": "random.Random seeded by first eight SHA256(case id) bytes; ACGT alphabet",
                "generation_edits_are_optimal_distances": False,
                "cases": [{key: value for key, value in case.items() if key not in {"x", "y"}}
                          for case in cases]}
    (directory / "cases.json").write_text(json.dumps(metadata, indent=2) + "\n")


def command_for(args, mode, bound):
    return [str(args.binary), "--mode", mode, "--bound", bound,
            "--width", str(args.width), "--block", str(args.block),
            "--seed-k", str(args.seed_k), "--epsilon", "0",
            "--order-block", str(args.order_block),
            "--order-max-cells", str(args.order_max_cells)]


def run_batch(args, setting, cases, repeat, batch_index, phase):
    # This is an in-memory adapter only; run_bound_bench.py stays unchanged.
    previous.command_for = command_for
    rows, batch, failures = previous.run_batch(args, setting, cases, repeat, batch_index, phase)
    for row in rows:
        label = f"{row['case_id']} {setting} {phase} {repeat}"
        for metric in ORDERED_METRICS:
            value = row.get(metric)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                failures.append(f"{label}: invalid ordered metric {metric}={value!r}")
        fallback = row.get("ordered_fallback")
        if not isinstance(fallback, str):
            failures.append(f"{label}: ordered_fallback must be a string reason (empty if none)")
        row["ordered_fallback_count"] = int(bool(fallback))
        row["combined_alignment_probe_cells"] = row.get("dp_cells", 0) + row.get("probe_cells", 0)
        row["combined_work_cells"] = row["combined_alignment_probe_cells"] + row.get("ordered_cells", 0)
    return rows, batch, failures


def scopes_for(cases):
    primary = [case for case in cases if case["suite"] != "fresh_v3"]
    scopes = [("all", cases), ("primary118", primary)]
    for suite in ("original", "heldout", "fresh_v3"):
        selected = [case for case in cases if case["suite"] == suite]
        scopes.append((suite, selected))
        for scenario in sorted({case["scenario"] for case in selected}):
            subset = [case for case in selected if case["scenario"] == scenario]
            scopes.append((f"{suite}/{scenario}", subset))
            if suite != "original":
                for gap in sorted({case["generation"]["indel_length"] for case in subset}):
                    scopes.append((f"{suite}/{scenario}/g{gap}",
                                   [case for case in subset if case["generation"]["indel_length"] == gap]))
    return scopes


def summarize(args, cases, records, batches, failures, elapsed_ms):
    timed = [row for row in records if row["phase"] == "timed"]
    metrics = previous.numeric_keys(timed)
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
    expected = len(cases) * len(SETTINGS) * args.repeats
    if len(timed) != expected:
        failures.append(f"expected {expected} timed outputs; got {len(timed)}")
    for case in cases:
        rows = [row for row in records if row["case_id"] == case["id"]]
        if len({row.get("distance") for row in rows}) != 1:
            failures.append(f"cross-setting/repetition distance mismatch for {case['id']}")
    aggregates, bound_comparisons, predictor_comparisons = [], [], []
    for scope, selected in scopes_for(cases):
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
                                                      "median_of_case_medians": statistics.median(values),
                                                      "max_of_case_medians": max(values)}
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
            result["case_change_counts"] = {}
            for metric in ("dp_cells", "rounds", "combined_work_cells", "total_ms"):
                counts = {"lower": 0, "equal": 0, "higher": 0}
                for case in selected:
                    litem = lookup.get((case["id"], *baseline))
                    ritem = lookup.get((case["id"], *candidate))
                    if litem is None or ritem is None or metric not in litem["metrics"] or metric not in ritem["metrics"]:
                        continue
                    lvalue = litem["metrics"][metric]["median"]
                    rvalue = ritem["metrics"][metric]["median"]
                    counts["lower" if rvalue < lvalue else "higher" if rvalue > lvalue else "equal"] += 1
                result["case_change_counts"][metric] = counts
            base_dp = result["metrics"].get("dp_ms", {}).get("baseline_sum_median", 0)
            cand_dp = result["metrics"].get("dp_ms", {}).get("candidate_sum_median", 0)
            result["observed_break_even"] = {
                "dp_time_saved_ms": base_dp - cand_dp,
                "added_other_time_ms": (b - cand_dp) - (a - base_dp),
                "net_total_time_saved_ms": a - b,
                "interpretation": "differences of sums of per-case median component times; preliminary shared-runtime observation, not a predictive performance model"}
            return result

        for mode in MODES:
            for baseline, candidate in (("length", "qgram"), ("length", "ordered"), ("qgram", "ordered")):
                bound_comparisons.append(compare((mode, baseline), (mode, candidate)))
        for bound in BOUNDS:
            for baseline, candidate in (("linear", "cost"), ("anchor-linear", "anchor-cost")):
                predictor_comparisons.append(compare((baseline, bound), (candidate, bound)))
    diagnostic = [row for row in timed if row["case_id"] == "compensating_indels-n4096-s17" and row["repeat"] == 0]
    fallback_reasons = {}
    for row in timed:
        if row["repeat"] == 0 and row.get("ordered_fallback"):
            reason = row["ordered_fallback"]
            fallback_reasons[reason] = fallback_reasons.get(reason, 0) + 1
    return {"schema_version": 3, "status": "fail" if failures else "pass",
            "configuration": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "fixed_configuration": {"modes": MODES, "bounds": BOUNDS, "epsilon": 0,
                "fresh_source_length": FRESH_LENGTH, "fresh_seeds": FRESH_SEEDS,
                "fresh_indel_lengths": FRESH_GAPS, "fresh_backgrounds": ["random", "repeated tract"],
                "warmups_per_setting": 1},
            "case_count": len(cases), "timed_result_count": len(timed),
            "result_count_including_warmup": len(records), "benchmark_wall_ms": elapsed_ms,
            "validation": {"failures": failures, "failure_count": len(failures),
                "evidence": "All CIGARs replay; exact certificates equal final scores; cross-setting/repetition agreement. Forty primary cases have an analytic optimum from a pure insertion/deletion construction. No independent full-DP oracle runs in this harness.",
                "analytic_reference_case_count": sum("analytic_optimal_distance" in case for case in cases),
                "trace_check": "Each finite trace lower/upper brackets the final certified distance; not independent proof of optimum."},
            "timing_policy": {
                "primary": "per-case median internal total_ms over five repetitions by default",
                "aggregation": "paired sums of per-case medians on identical cases, separately by suite/scenario",
                "order": "12 mode/bound settings; rotate by 5*repeat positions, reverse odd repetitions",
                "warmup": "one complete untimed batch per setting, retained in raw outputs and excluded from medians",
                "clock": "engine elapsed wall-clock milliseconds, not CPU cycles",
                "process": "fresh process per batch; wall_ms separately includes process startup and I/O",
                "overheads": "all engine preprocessing, prediction and certification included in total_ms",
                "external_baselines": "none; no external-library or state-of-art speed claim",
                "combined_work_cells": "main DP plus predictor probe plus ordered-index table cells (including its terminal row); unequal per-cell cost and exclusion of qgram/anchor/certificate operations make this a partial work proxy, not elapsed-time equivalence",
                "ordered_fallback_count": "number of per-case medians with a nonempty ordered_fallback reason; warmups and repeated runs not multiply counted",
                "post_first_optimal_dp_cells": "cells in rounds after the first upper equal to the final exact score, identified retrospectively; excludes the first successful round"},
            "metric_names": metrics, "batches": batches, "per_case": per_case,
            "diagnostic_D818_case": diagnostic,
            "ordered_fallback_reasons_first_timed_repetition": fallback_reasons,
            "aggregates": aggregates, "paired_bound_comparisons": bound_comparisons,
            "paired_predictor_comparisons": predictor_comparisons}


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "build" / "certalign")
    parser.add_argument("--primary-dir", type=Path, default=root / "results" / "bound_benchmark")
    parser.add_argument("--original-dir", type=Path, default=root / "results" / "benchmark")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "ordered_benchmark")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--width", type=int, default=8)
    parser.add_argument("--block", type=int, default=64)
    parser.add_argument("--seed-k", type=int, default=11)
    parser.add_argument("--order-block", type=int, default=8)
    parser.add_argument("--order-max-cells", type=int, default=4_000_000)
    parser.add_argument("--timeout", type=float, default=180.)
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args(argv)
    if min(args.repeats, args.width, args.block, args.seed_k, args.order_block, args.order_max_cells) < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("counts, sizes and timeout must be positive")
    for key in ("binary", "primary_dir", "original_dir", "output_dir"):
        setattr(args, key, getattr(args, key).resolve())
    if args.output_dir in (args.primary_dir, args.original_dir):
        parser.error("output-dir must differ from input directories")
    if not args.generate_only and not args.binary.is_file():
        parser.error(f"binary does not exist: {args.binary}")
    predictor_sha = previous.sha256(root / "include" / "predictor.hpp")
    if predictor_sha != previous.PREDICTOR_SHA256:
        parser.error("predictor.hpp differs from the fixed v1 predictor")
    cases = load_primary(args.primary_dir, args.original_dir)
    cases += [make_fresh(seed, gap, repeat) for seed in FRESH_SEEDS for gap in FRESH_GAPS for repeat in (False, True)]
    if len(cases) != 148 or len({case["id"] for case in cases}) != 148:
        raise ValueError("invalid fixed v3 suite")
    write_inputs(args.output_dir, cases)
    fingerprints = {"predictor_sha256": predictor_sha, "primary_118_pairs_sha256": PRIMARY_SHA256,
                    "original_pairs_sha256": previous.ORIGINAL_PAIRS_SHA256,
                    "pairs_tsv_sha256": previous.sha256(args.output_dir / "pairs.tsv")}
    if args.generate_only:
        print(json.dumps({"case_count": 148, "primary": 118, "fresh": 30, **fingerprints,
                          "output_dir": str(args.output_dir)}, indent=2))
        return 0
    binary_sha = previous.sha256(args.binary)
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
    summary.update(engine_binary_sha256=binary_sha, **fingerprints)
    previous.save_results(args.output_dir, records, summary)
    print(json.dumps({"status": summary["status"], "case_count": len(cases),
                      "timed_result_count": summary["timed_result_count"], "wall_seconds": elapsed_ms / 1000,
                      "failure_count": len(failures), "output_dir": str(args.output_dir)}, indent=2))
    for failure in failures[:10]:
        print(f"FAIL: {failure}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
