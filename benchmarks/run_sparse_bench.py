#!/usr/bin/env python3
"""Frozen v4 sparse-certificate ablation, preserving all 148 prior pairs.

Four predictors by four suffix bounds, five timed repetitions and one warmup.
The new thirty compensated-indel pairs are fixed before measuring v4 outputs.
The primary milestone compares each sparse variant with plain linear/length.
No file in an earlier benchmark harness is modified by this harness.
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
BOUNDS = ("length", "qgram", "ordered", "sparse")
SETTINGS = tuple((mode, bound) for mode in MODES for bound in BOUNDS)
PRIMARY_SHA256 = "80e6c6e100848eb13a364a9395fec7e82191cbbe69ca8d47444c619b78249579"
FRESH_SEEDS = (401, 409, 419, 421, 431)
FRESH_GAPS = (47, 193, 511)
FRESH_LENGTH = 3073
SUITES = ("original63", "v2extension55", "v3fresh30", "NEWfreshv4_30")
REQUIRED_EXTRA_METRICS = ("bound_preprocess_ms", "bound_bytes", "ordered_cells", "ordered_bytes", "ordered_preprocess_ms", "ordered_calls", "sparse_queries", "sparse_expansions", "sparse_cache_hits", "sparse_guard_rejections", "sparse_cached_states", "sparse_tiles", "sparse_bytes", "sparse_peak_stack", "sparse_preprocess_ms")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def load_primary(directory, original_directory):
    original = previous.load_original(original_directory)
    pairs_path = directory / "pairs.tsv"
    if previous.sha256(pairs_path) != PRIMARY_SHA256:
        raise ValueError("primary 148 pairs differ from published v3 inputs")
    metadata = json.loads((directory / "cases.json").read_text())["cases"]
    by_id = {case["id"]: case for case in metadata}
    cases = []
    renamed = {"original": "original63", "heldout": "v2extension55", "fresh_v3": "v3fresh30"}
    for line in pairs_path.read_text().splitlines():
        case_id, x, y = line.split("\t")
        case = dict(by_id[case_id], x=x, y=y)
        if hashlib.sha256((x + "\n" + y).encode()).hexdigest() != case["input_sha256"]:
            raise ValueError(f"metadata digest mismatch for {case_id}")
        if (len(x), len(y)) != (case["x_length"], case["y_length"]):
            raise ValueError(f"metadata length mismatch for {case_id}")
        case["prior_suite"] = case["suite"]
        case["suite"] = renamed[case["suite"]]
        case["evidence_status"] = "known from prior experiments; not a fresh test"
        cases.append(case)
    if len(cases) != 148 or len(by_id) != 148:
        raise ValueError("primary suite must contain exactly 148 cases")
    existing = {case["id"]: (case["x"], case["y"]) for case in cases}
    if any(existing.get(case["id"]) != (case["x"], case["y"]) for case in original):
        raise ValueError("primary suite does not preserve all 63 original pairs")
    if sum("analytic_optimal_distance" in case for case in cases) != 40:
        raise ValueError("primary suite must preserve forty analytic references")
    return cases


def make_fresh(seed, gap, repeat_background):
    scenario = "repeat_compensating" if repeat_background else "compensating"
    case_id = f"fresh-v4-{scenario}-n{FRESH_LENGTH}-g{gap}-s{seed}"
    derived_seed = int.from_bytes(hashlib.sha256(case_id.encode()).digest()[:8], "big")
    rng = random.Random(derived_seed)
    x = previous.dna(rng, FRESH_LENGTH)
    details = {"indel_length": gap, "inserted_bases": gap, "deleted_bases": gap,
               "substitutions": 0, "repeat_background": repeat_background,
               "background": "repeated" if repeat_background else "random"}
    if repeat_background:
        motif = previous.dna(rng, 4)
        start, stop = FRESH_LENGTH // 4, FRESH_LENGTH * 3 // 4
        tract = (motif * ((stop - start + 3) // 4))[:stop - start]
        x = x[:start] + tract + x[stop:]
        details.update(motif=motif, repeat_interval=[start, stop])
    deletion_position = FRESH_LENGTH // 4 + rng.randrange(64)
    short = x[:deletion_position] + x[deletion_position + gap:]
    insertion_position = len(short) * 3 // 4 + rng.randrange(64)
    inserted = previous.dna(rng, gap)
    y = short[:insertion_position] + inserted + short[insertion_position:]
    details.update(deletion_position=deletion_position,
                   insertion_position_after_deletion=insertion_position)
    return {"id": case_id, "x": x, "y": y, "suite": "NEWfreshv4_30", "scenario": scenario,
            "evidence_status": "new v4 inputs fixed before v4 timing; no tuning on outcomes",
            "source_length": FRESH_LENGTH, "seed": seed, "derived_seed": derived_seed,
            "x_length": len(x), "y_length": len(y), "generation": details,
            "input_sha256": hashlib.sha256((x + "\n" + y).encode()).hexdigest()}


def write_frozen(path, contents):
    if path.exists() and path.read_text() != contents:
        raise ValueError(f"refusing to change previously frozen input/configuration: {path}")
    path.write_text(contents)


def freeze_inputs(args, cases, predictor_sha):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_frozen(args.output_dir / "pairs.tsv", "".join(
        f"{case['id']}\t{case['x']}\t{case['y']}\n" for case in cases))
    metadata = {"schema_version": 4,
                "primary_policy": "published 148 v3 pairs loaded unchanged and SHA256 verified; all are now known",
                "original_policy": "original 63 v1 pairs additionally verified by published SHA256",
                "fresh_policy": "30 new compensated-indel pairs; fixed seeds, gaps, backgrounds before v4 timing",
                "fresh_generator": "random.Random seeded by first eight SHA256(case id) bytes; ACGT alphabet",
                "generation_edits_are_optimal_distances": False,
                "cases": [{key: value for key, value in case.items() if key not in {"x", "y"}} for case in cases]}
    write_frozen(args.output_dir / "cases.json", json.dumps(metadata, indent=2) + "\n")
    config = {"schema_version": 4, "modes": MODES, "bounds": BOUNDS,
              "width": args.width, "block": args.block, "seed_k": args.seed_k, "epsilon": 0,
              "order_block": args.order_block, "order_max_cells": args.order_max_cells,
              "sparse_max_work": args.sparse_max_work, "bound_memory_bytes": args.bound_memory,
              "repetitions": args.repeats, "warmups_per_setting": 1,
              "setting_order": "rotate by five positions per repetition; reverse odd repetitions",
              "fresh_source_length": FRESH_LENGTH, "fresh_seeds": FRESH_SEEDS,
              "fresh_indel_lengths": FRESH_GAPS, "fresh_backgrounds": ["random", "repeated"],
              "primary_148_pairs_sha256": PRIMARY_SHA256,
              "pairs_tsv_sha256": previous.sha256(args.output_dir / "pairs.tsv"),
              "predictor_sha256": predictor_sha,
              "primary_milestone": "each sparse variant versus plain linear predictor and length bound on matched cases; fastest sparse winner selected at aggregate scope, not per case",
              "analytic_reference_count": 40}
    config["configuration_sha256"] = digest_json(config)
    write_frozen(args.output_dir / "frozen_configuration.json", json.dumps(config, indent=2) + "\n")
    return config


def command_for(args, mode, bound):
    return [str(args.binary), "--mode", mode, "--bound", bound,
            "--width", str(args.width), "--block", str(args.block),
            "--seed-k", str(args.seed_k), "--epsilon", "0",
            "--order-block", str(args.order_block), "--order-max-cells", str(args.order_max_cells),
            "--sparse-max-work", str(args.sparse_max_work), "--bound-memory", str(args.bound_memory)]


def run_batch(args, setting, cases, repeat, batch_index, phase):
    previous.command_for = command_for  # in-memory adapter; never modify earlier harness files
    rows, batch, failures = previous.run_batch(args, setting, cases, repeat, batch_index, phase)
    for row in rows:
        label = f"{row['case_id']} {setting} {phase} {repeat}"
        for metric in REQUIRED_EXTRA_METRICS:
            value = row.get(metric)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                failures.append(f"{label}: invalid metric {metric}={value!r}")
        for metric, value in row.items():
            if metric.startswith("sparse_") and isinstance(value, (float, int)) and not isinstance(value, bool):
                if not math.isfinite(value) or value < 0:
                    failures.append(f"{label}: invalid sparse metric {metric}={value!r}")
        for name in ("bound_fallback", "ordered_fallback", "sparse_fallback"):
            value = row.get(name, "")
            if not isinstance(value, str):
                failures.append(f"{label}: {name} must be a string")
            row[name + "_count"] = int(bool(value))
        row["combined_alignment_probe_cells"] = row.get("dp_cells", 0) + row.get("probe_cells", 0)
        row["combined_work_units"] = row["combined_alignment_probe_cells"] + row.get("ordered_cells", 0) + row.get("sparse_expansions", 0)
    return rows, batch, failures


def scopes_for(cases):
    scopes = [("all178", cases), ("primary148", [case for case in cases if case["suite"] != "NEWfreshv4_30"])]
    for suite in SUITES:
        selected = [case for case in cases if case["suite"] == suite]
        scopes.append((suite, selected))
        for scenario in sorted({case["scenario"] for case in selected}):
            subset = [case for case in selected if case["scenario"] == scenario]
            scopes.append((f"{suite}/{scenario}", subset))
            if suite != "original63":
                for gap in sorted({case["generation"]["indel_length"] for case in subset}):
                    scopes.append((f"{suite}/{scenario}/g{gap}", [case for case in subset if case["generation"]["indel_length"] == gap]))
        if suite in ("v3fresh30", "NEWfreshv4_30"):
            for gap in sorted({case["generation"]["indel_length"] for case in selected}):
                scopes.append((f"{suite}/all_backgrounds/g{gap}", [case for case in selected if case["generation"]["indel_length"] == gap]))
    return scopes


def summarize(args, cases, records, batches, failures, elapsed_ms, frozen):
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
        item["fallback_reasons"] = {key: sorted({row.get(key, "") for row in rows} - {""})
                                    for key in ("bound_fallback", "ordered_fallback", "sparse_fallback")}
        if len(rows) != args.repeats:
            failures.append(f"{case_id} {mode}/{bound}: expected {args.repeats} repetitions; got {len(rows)}")
        lookup[(case_id, mode, bound)] = item
        per_case.append(item)
    expected = len(cases) * len(SETTINGS) * args.repeats
    if len(timed) != expected:
        failures.append(f"expected {expected} timed outputs; got {len(timed)}")
    for case in cases:
        rows = [row for row in records if row["case_id"] == case["id"]]
        if len({row.get("distance") for row in rows}) != 1:
            failures.append(f"cross-setting/repetition distance mismatch for {case['id']}")
    aggregates, comparisons, milestones, sparse_diagnostics = [], [], [], []
    scopes = scopes_for(cases)
    for scope, selected in scopes:
        sums = {}
        for mode, bound in SETTINGS:
            items = [lookup[(case["id"], mode, bound)] for case in selected if (case["id"], mode, bound) in lookup]
            aggregated = {"scope": scope, "mode": mode, "bound": bound, "case_count": len(items), "metrics": {}}
            for metric in metrics:
                values = [item["metrics"][metric]["median"] for item in items if metric in item["metrics"]]
                if values:
                    aggregated["metrics"][metric] = {"sum_of_case_medians": sum(values),
                                                      "median_of_case_medians": statistics.median(values),
                                                      "max_of_case_medians": max(values)}
            total = aggregated["metrics"].get("dp_cells", {}).get("sum_of_case_medians", 0)
            post = aggregated["metrics"].get("post_first_optimal_dp_cells", {}).get("sum_of_case_medians", 0)
            aggregated["post_first_optimal_fraction_of_total_cells"] = post / total if total else None
            sums[(mode, bound)] = aggregated
            aggregates.append(aggregated)

        def compare(baseline, candidate, kind):
            left, right = sums[baseline], sums[candidate]
            result = {"scope": scope, "kind": kind, "baseline_mode": baseline[0], "baseline_bound": baseline[1],
                      "candidate_mode": candidate[0], "candidate_bound": candidate[1], "matched_case_count": len(selected), "metrics": {}}
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
            for metric in ("dp_cells", "rounds", "combined_work_units", "total_ms"):
                counts = {"lower": 0, "equal": 0, "higher": 0}
                for case in selected:
                    litem, ritem = lookup.get((case["id"], *baseline)), lookup.get((case["id"], *candidate))
                    if litem is None or ritem is None or metric not in litem["metrics"] or metric not in ritem["metrics"]:
                        continue
                    lv, rv = litem["metrics"][metric]["median"], ritem["metrics"][metric]["median"]
                    counts["lower" if rv < lv else "higher" if rv > lv else "equal"] += 1
                result["case_change_counts"][metric] = counts
            return result

        for mode in MODES:
            for base, candidate in (("length", "qgram"), ("length", "ordered"), ("qgram", "ordered"),
                                    ("ordered", "sparse"), ("qgram", "sparse"), ("length", "sparse")):
                comparisons.append(compare((mode, base), (mode, candidate), "same_predictor_bound_ablation"))
            milestone = compare(("linear", "length"), (mode, "sparse"), "primary_milestone")
            comparisons.append(milestone)
            milestones.append(milestone)
        for bound in BOUNDS:
            for base, candidate in (("linear", "cost"), ("anchor-linear", "anchor-cost")):
                comparisons.append(compare((base, bound), (candidate, bound), "same_bound_predictor_ablation"))
        # Winner uses one setting across a scope, never an impossible best-setting-per-pair oracle.
        available = [(mode, "sparse") for mode in MODES if "total_ms" in sums[(mode, "sparse")]["metrics"]]
        if available:
            winner = min(available, key=lambda s: sums[s]["metrics"]["total_ms"]["sum_of_case_medians"])
            best = compare(("linear", "length"), winner, "fastest_sparse_milestone")
            best["selection_caveat"] = "four candidate settings selected using this scope's timings; descriptive, not an independently validated model-selection result"
            milestones.append(best)
        for mode in MODES:
            counts = {"same_dp_cells_as_dense": 0, "more_dp_cells_than_dense": 0, "fewer_dp_cells_than_dense": 0,
                      "more_dp_cells_with_sparse_fallback": 0, "sparse_fallback_cases": 0, "dense_fallback_cases": 0}
            fallback_reasons = {}
            for case in selected:
                sparse, dense = lookup.get((case["id"], mode, "sparse")), lookup.get((case["id"], mode, "ordered"))
                if sparse is None or dense is None:
                    continue
                sdp, ddp = sparse["metrics"]["dp_cells"]["median"], dense["metrics"]["dp_cells"]["median"]
                counts["same_dp_cells_as_dense" if sdp == ddp else "more_dp_cells_than_dense" if sdp > ddp else "fewer_dp_cells_than_dense"] += 1
                reasons = sparse["fallback_reasons"]["sparse_fallback"]
                if reasons:
                    counts["sparse_fallback_cases"] += 1
                    counts["more_dp_cells_with_sparse_fallback"] += int(sdp > ddp)
                    for reason in reasons:
                        fallback_reasons[reason] = fallback_reasons.get(reason, 0) + 1
                counts["dense_fallback_cases"] += int(bool(dense["fallback_reasons"]["ordered_fallback"]))
            comparison = compare((mode, "ordered"), (mode, "sparse"), "sparse_diagnostic")
            sparse_diagnostics.append({"scope": scope, "mode": mode, "case_count": len(selected), **counts,
                                       "sparse_fallback_reasons": fallback_reasons,
                                       "preprocessing_metrics": {key: val for key, val in comparison["metrics"].items() if "preprocess" in key or "bytes" in key or key.startswith("sparse_")},
                                       "causality_note": "extra DP work with fallback is an observed association; cap causality is not established without an uncapped rerun"})
    diagnostic = [row for row in timed if row["case_id"] == "compensating_indels-n4096-s17" and row["repeat"] == 0]
    return {"schema_version": 4, "status": "fail" if failures else "pass", "frozen_configuration": frozen,
            "configuration": {key: str(value) if isinstance(value, Path) else value for key, value in vars(args).items()},
            "case_count": len(cases), "timed_result_count": len(timed), "result_count_including_warmup": len(records),
            "benchmark_wall_ms": elapsed_ms,
            "validation": {"failure_count": len(failures), "failures": failures,
                "evidence": "Every CIGAR replays; all exact certificates equal scores; scores agree across settings and repetitions; forty old cases have independent analytic optima. No full-DP oracle runs in this benchmark harness.",
                "analytic_reference_case_count": sum("analytic_optimal_distance" in case for case in cases),
                "trace_check": "every finite trace interval contains final certified distance; not independent proof of optimum"},
            "timing_policy": {"primary": "sum of per-case median internal total_ms over five repetitions by default, on identical cases",
                "order": "16 settings rotate by five positions per repetition and reverse odd repetitions",
                "warmup": "one complete batch per setting, excluded from medians but retained in raw outputs",
                "overheads": "total_ms includes all engine prediction, preprocessing, lazy expansion, DP, certificate and traceback overhead",
                "wall_ms": "batch wall includes process startup and I/O; separately reported",
                "bound_memory": "64 MiB default shared bound-memory cap; explicitly passed to all settings",
                "work_proxy": "combined_work_units = main DP cells + predictor probe cells + ordered table cells + sparse state expansions; costs differ and this excludes other work",
                "fallback_counts": "number of cases with fallback, not multiplied by repetitions or warmups",
                "post_first_optimal": "DP cells in rounds after first upper equal to final exact score; retrospective diagnostic",
                "limits": "synthetic cases and one shared runtime; no external library or state-of-art comparison; previous 148 pairs are known data"},
            "metric_names": metrics, "batches": batches, "per_case": per_case, "aggregates": aggregates,
            "paired_comparisons": comparisons, "primary_milestones": milestones, "sparse_diagnostics": sparse_diagnostics,
            "diagnostic_D818_case": diagnostic}


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "build" / "certalign")
    parser.add_argument("--primary-dir", type=Path, default=root / "results" / "ordered_benchmark")
    parser.add_argument("--original-dir", type=Path, default=root / "results" / "benchmark")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "sparse_benchmark")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--width", type=int, default=8)
    parser.add_argument("--block", type=int, default=64)
    parser.add_argument("--seed-k", type=int, default=11)
    parser.add_argument("--order-block", type=int, default=8)
    parser.add_argument("--order-max-cells", type=int, default=4_000_000)
    parser.add_argument("--sparse-max-work", type=int, default=4_000_000)
    parser.add_argument("--bound-memory", type=int, default=64 * 1024 * 1024)
    parser.add_argument("--timeout", type=float, default=600.)
    parser.add_argument("--generate-only", action="store_true")
    args = parser.parse_args(argv)
    if min(args.repeats, args.width, args.block, args.seed_k, args.order_block, args.order_max_cells,
           args.sparse_max_work, args.bound_memory) < 1 or not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("counts, sizes and timeout must be positive")
    for key in ("binary", "primary_dir", "original_dir", "output_dir"):
        setattr(args, key, getattr(args, key).resolve())
    if args.output_dir in (args.primary_dir, args.original_dir):
        parser.error("output directory must differ from input directories")
    predictor_sha = previous.sha256(root / "include" / "predictor.hpp")
    if predictor_sha != previous.PREDICTOR_SHA256:
        parser.error("predictor.hpp differs from the fixed v1 predictor")
    cases = load_primary(args.primary_dir, args.original_dir)
    cases += [make_fresh(seed, gap, repeat) for seed in FRESH_SEEDS for gap in FRESH_GAPS for repeat in (False, True)]
    if len(cases) != 178 or len({case["id"] for case in cases}) != 178:
        raise ValueError("invalid fixed v4 suite")
    frozen = freeze_inputs(args, cases, predictor_sha)
    print(json.dumps({"phase": "inputs_frozen", "case_count": 178, "primary": 148, "NEWfreshv4": 30,
                      "configuration_sha256": frozen["configuration_sha256"],
                      "pairs_tsv_sha256": frozen["pairs_tsv_sha256"], "output_dir": str(args.output_dir)}, indent=2), flush=True)
    if args.generate_only:
        return 0
    if not args.binary.is_file():
        parser.error(f"binary does not exist: {args.binary}")
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
    if previous.sha256(args.binary) != binary_sha:
        failures.append("engine binary changed during timings")
    if previous.sha256(root / "include" / "predictor.hpp") != predictor_sha:
        failures.append("predictor source changed during timings")
    summary = summarize(args, cases, records, batches, failures, elapsed_ms, frozen)
    summary.update(engine_binary_sha256=binary_sha, predictor_sha256=predictor_sha,
                   primary_148_pairs_sha256=PRIMARY_SHA256, pairs_tsv_sha256=frozen["pairs_tsv_sha256"])
    previous.save_results(args.output_dir, records, summary)
    print(json.dumps({"status": summary["status"], "case_count": len(cases),
                      "timed_result_count": summary["timed_result_count"], "wall_seconds": elapsed_ms / 1000,
                      "failure_count": len(failures), "output_dir": str(args.output_dir)}, indent=2))
    for failure in failures[:10]:
        print(f"FAIL: {failure}", file=sys.stderr)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
