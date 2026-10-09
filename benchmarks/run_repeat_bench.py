#!/usr/bin/env python3
"""Frozen v5 repeat validation of one previously selected v4 setting.

Generate/freeze inputs first; timing requires independently computed oracle scores.
The primary endpoint is the 60-case replication suite. Controlled repeats and
matched random controls are secondary, never pooled into that endpoint.
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
from collections import Counter
from pathlib import Path

import run_bench as original
import run_bound_bench as previous
import run_sparse_bench as sparse_v4

BOUNDS = ("length", "qgram", "ordered", "sparse")
SETTINGS = tuple(("linear", bound) for bound in BOUNDS)
EXPECTED_BINARY_SHA256 = "4a492041691bd03ab0577152cfe0aa056a19f8f4be18551e43538d4a8d476a25"
REPLICATION_SEEDS = tuple(range(501, 521))
LOCUS_SEEDS = (601, 607)
PERIODS = (2, 3, 4, 5, 6)
SUITES = ("fresh_replication60", "controlled_repeats160", "matched_random160")
BOOTSTRAP_SEED = 20261009
BOOTSTRAP_REPLICATES = 9999
CASE_ORDER_SEED = 2026100901
ALPHABET = "ACGT"


def seed_for(label):
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def minimal_period(motif):
    return next(p for p in range(1, len(motif) + 1)
                if all(motif[i] == motif[i - p] for i in range(p, len(motif))))


def primitive_repeat_period(motif):
    return next(p for p in range(1, len(motif) + 1)
                if len(motif) % p == 0 and motif == motif[:p] * (len(motif) // p))


def substitution_mask(rng, start, stop, probability):
    """Bernoulli mask with a nonzero cyclic alphabet offset at each position."""
    return [[i, rng.randrange(1, 4)] for i in range(start, stop)
            if rng.random() < probability]


def apply_mask(sequence, mask):
    values = list(sequence)
    for position, offset in mask:
        values[position] = ALPHABET[(ALPHABET.index(values[position]) + offset) % 4]
    return "".join(values)


def make_replication():
    cases = []
    for n in (256, 1024, 4096):
        for seed in REPLICATION_SEEDS:
            case = original.generate_case(n, seed, "tandem_repeats")
            case.update(suite=SUITES[0], cluster=f"rep-seed-{seed}",
                        evidence_status="fresh seed replication of the original tandem_repeats generator",
                        generator="unmodified run_bench.generate_case")
            motif = case["generation"]["motif"]
            case["generation"].update(raw_motif_length=len(motif),
                minimal_string_period=minimal_period(motif),
                effective_primitive_period=primitive_repeat_period(motif))
            cases.append(case)
    return cases


def make_locus(n, q, seed):
    """Create both backgrounds jointly so all edit/noise masks are matched."""
    label = f"v5-locus-n{n}-q{q}-s{seed}"
    rng = random.Random(seed_for(label))
    random_source = original.dna(rng, n)
    motif = original.dna(rng, q)
    while primitive_repeat_period(motif) != q:
        motif = original.dna(rng, q)
    copies = (n // 2) // q
    tract_length = copies * q
    start = (n - tract_length) // 2
    stop = start + tract_length
    repeat_source = random_source[:start] + motif * copies + random_source[stop:]
    # Both backgrounds share the same flanks and exact noise coordinates/offsets.
    source_mask = substitution_mask(random.Random(seed_for(label + "/source-noise")), start, stop, .01)
    position = start + (copies // 2) * q
    cases = []
    for percent in (2, 8):
        delta_copies = max(1, (n * percent) // (100 * q))
        gap = delta_copies * q
        inserted_random = original.dna(random.Random(seed_for(f"{label}/insert/{percent}")), gap)
        for direction in ("expand", "contract"):
            target_length = n + gap if direction == "expand" else n - gap
            shared_target_mask = substitution_mask(
                random.Random(seed_for(f"{label}/target-noise/{percent}/{direction}")),
                0, target_length, .004)
            for regime in ("perfect", "imperfect"):
                active_source_mask = source_mask if regime == "imperfect" else []
                active_target_mask = shared_target_mask if regime == "imperfect" else []
                pair_id = f"{label}-d{percent}-{direction}-{regime}"
                for background, base_source, suite in (
                        ("repeat", repeat_source, SUITES[1]),
                        ("random", random_source, SUITES[2])):
                    x = apply_mask(base_source, active_source_mask)
                    inserted = motif * delta_copies if background == "repeat" else inserted_random
                    if direction == "expand":
                        y_pre_noise = x[:position] + inserted + x[position:]
                    else:
                        # Delete the actual source, preserving inherited interruptions elsewhere.
                        y_pre_noise = x[:position] + x[position + gap:]
                    y = apply_mask(y_pre_noise, active_target_mask)
                    details = {
                        "background": background, "motif": motif, "primitive_period": q,
                        "motif_period_definition": "smallest divisor p of motif length such that motif is an exact power of its length-p prefix",
                        "minimal_string_period": minimal_period(motif),
                        "tract_interval": [start, stop], "tract_copies": copies,
                        "position": position, "delta_percent": percent,
                        "indel_length": gap, "delta_copies": delta_copies,
                        "direction": direction, "regime": regime,
                        "inserted_bases": gap if direction == "expand" else 0,
                        "deleted_bases": gap if direction == "contract" else 0,
                        "source_interruption_probability": .01 if regime == "imperfect" else 0,
                        "source_interruption_mask": active_source_mask,
                        "source_interruption_count": len(active_source_mask),
                        "removed_source_interruptions": sum(position <= i < position + gap
                            for i, _ in active_source_mask) if direction == "contract" else 0,
                        "target_substitution_probability": .004 if regime == "imperfect" else 0,
                        "target_substitution_mask": active_target_mask,
                        "target_substitution_count": len(active_target_mask),
                        "substitution_offset_rule": "replace ACGT index b by (b+offset) mod4, offset in1..3",
                        "inserted_sequence": inserted if direction == "expand" else "",
                        "deleted_sequence": x[position:position + gap] if direction == "contract" else "",
                        "pre_noise_target_sha256": hashlib.sha256(y_pre_noise.encode()).hexdigest(),
                        "regime_caveat": "imperfect jointly introduces source interruptions and target noise; effects not separated",
                    }
                    case = {"id": f"{pair_id}-{background}", "x": x, "y": y,
                            "suite": suite, "scenario": "repeat_copy_change" if background == "repeat" else "matched_random_indel",
                            "source_length": n, "seed": seed, "derived_seed": seed_for(label),
                            "x_length": len(x), "y_length": len(y), "cluster": label,
                            "matched_pair_id": pair_id, "generation": details,
                            "generation_sha256": digest_json(details),
                            "evidence_status": "fresh controlled transfer set fixed before v5 timing",
                            "input_sha256": hashlib.sha256((x + "\n" + y).encode()).hexdigest()}
                    if regime == "perfect":
                        case["analytic_optimal_distance"] = gap
                        case["analytic_optimum_reason"] = "length difference equals a constructed single insertion/deletion path"
                    cases.append(case)
    return cases


def make_cases():
    cases = make_replication()
    for n in (1024, 4096):
        for q in PERIODS:
            for seed in LOCUS_SEEDS:
                cases.extend(make_locus(n, q, seed))
    if Counter(case["suite"] for case in cases) != dict(zip(SUITES, (60, 160, 160))):
        raise ValueError("unexpected frozen suite counts")
    if len({case["id"] for case in cases}) != 380:
        raise ValueError("duplicate case ids")
    if sum("analytic_optimal_distance" in case for case in cases) != 160:
        raise ValueError("expected160 analytic pure-indel references")
    return cases


def source_hashes(root):
    paths = sorted([*root.glob("include/*.hpp"), *root.glob("src/*.cpp"),
                    root / "REPEAT_VALIDATION_PLAN.md",
                    *(root / "benchmarks" / name for name in
                      ("run_bench.py", "run_bound_bench.py", "run_sparse_bench.py", "run_repeat_bench.py"))])
    return {str(path.relative_to(root)): previous.sha256(path) for path in paths}


def freeze_inputs(args, cases, root):
    args.output_dir.mkdir(parents=True, exist_ok=True)
    old_cases = json.loads((root / "results" / "sparse_benchmark" / "cases.json").read_text())["cases"]
    overlap = set(c["input_sha256"] for c in cases) & set(c["input_sha256"] for c in old_cases)
    if overlap:
        raise ValueError("fresh cases duplicate prior v4 input pairs")
    hashes = {}
    for case in cases:
        hashes.setdefault(case["input_sha256"], []).append(case["id"])
    duplicate_groups = [ids for ids in hashes.values() if len(ids) > 1]
    sparse_v4.write_frozen(args.output_dir / "pairs.tsv", "".join(
        f"{c['id']}\t{c['x']}\t{c['y']}\n" for c in cases))
    metadata = {"schema_version": 5, "generation_edits_are_optimal_distances": False,
                "suite_counts": dict(Counter(c["suite"] for c in cases)),
                "policy": "fixed fresh seeds and paired controls; engine and parameters unchanged from v4",
                "replication": "unmodified old generator, preserving its exact x/y and id on new seeds501..520",
                "controlled": "20 loci,8 variants each in both backgrounds; masks and geometry paired",
                "freshness_check": {"prior_case_count": len(old_cases), "prior_input_overlap_count": len(overlap),
                                    "fresh_exact_duplicate_groups": duplicate_groups},
                "coordinate_convention": "zero-based positions and half-open tract/deletion intervals",
                "cases": [{k: v for k, v in c.items() if k not in {"x", "y"}} for c in cases]}
    sparse_v4.write_frozen(args.output_dir / "cases.json", json.dumps(metadata, indent=2) + "\n")
    config = {
        "schema_version": 5, "candidate": "linear/sparse", "co_primary_baselines": ["linear/length", "linear/qgram"],
        "secondary_baseline": "linear/ordered", "settings": SETTINGS,
        "parameters": {"width": args.width, "block": args.block, "seed_k": args.seed_k, "epsilon": 0,
                       "order_block": args.order_block, "order_max_cells": args.order_max_cells,
                       "sparse_max_work": args.sparse_max_work, "bound_memory": args.bound_memory,
                       "main_dp_max_cells": 50_000_000},
        "repetitions": args.repeats, "warmups_per_setting": 1,
        "setting_order": "rotate by repetition mod4; reverse odd repetitions",
        "case_order": "random.Random(case_order_seed+repeat+1).shuffle; repeat=-1 for warmup; identical across settings",
        "case_order_seed": CASE_ORDER_SEED,
        "replication_seeds": REPLICATION_SEEDS, "replication_lengths": [256, 1024, 4096],
        "controlled_lengths": [1024, 4096], "primitive_periods": PERIODS, "locus_seeds": LOCUS_SEEDS,
        "delta_percent": [2, 8], "directions": ["expand", "contract"], "regimes": ["perfect", "imperfect"],
        "delta_rule": "max(1,floor(n*percent/(100*period)))*period",
        "primary_scope": SUITES[0],
        "primary_gate": "candidate sum of per-case median total_ms <=0.90 times EACH co-primary baseline AND each paired cluster-bootstrap95% interval upper ratio <1",
        "bootstrap": {"replicates": BOOTSTRAP_REPLICATES, "seed": BOOTSTRAP_SEED,
                      "method": "percentile95% with linear-interpolated2.5/97.5 percentiles; paired resampling of whole clusters",
                      "replication_cluster": "seed; each contains3 lengths", "controlled_cluster": "n/period/seed locus; includes all8 variants",
                      "scope_policy": "each suite separate; subgroups descriptive and unadjusted; no pooled endpoint",
                      "uncertainty_limit": "clusters of synthetic inputs; timings from one runtime; not cross-machine uncertainty"},
        "engine_binary_sha256": args.binary_sha256,
        "original_v4_engine_binary_sha256": EXPECTED_BINARY_SHA256,
        "portable_rebuilt_binary": args.allow_rebuilt_binary,
        "predictor_sha256": previous.PREDICTOR_SHA256,
        "source_sha256": source_hashes(root),
        "pairs_tsv_sha256": previous.sha256(args.output_dir / "pairs.tsv"),
        "metadata_sha256": previous.sha256(args.output_dir / "cases.json"),
        "analytic_reference_count": 160,
        "correctness": "independent full-DP oracle required for all380 pairs, plus all output CIGARs and round bounds",
        "selection_policy": "single linear/sparse candidate chosen on old results; no model selection or tuning on these380 outcomes",
    }
    config["configuration_sha256"] = digest_json(config)
    sparse_v4.write_frozen(args.output_dir / "frozen_configuration.json", json.dumps(config, indent=2) + "\n")
    return config


def load_oracle(path, frozen, cases):
    data = json.loads(path.read_text())
    if data.get("pairs_sha256") != frozen["pairs_tsv_sha256"]:
        raise ValueError("oracle pairs hash mismatch")
    if data.get("metadata_sha256") != frozen["metadata_sha256"]:
        raise ValueError("oracle metadata hash mismatch")
    scores = data.get("scores", {})
    if set(scores) != {c["id"] for c in cases}:
        raise ValueError("oracle score ids mismatch")
    for c in cases:
        d = scores[c["id"]]
        if isinstance(d, bool) or not isinstance(d, int) or not abs(c["x_length"] - c["y_length"]) <= d <= max(c["x_length"], c["y_length"]):
            raise ValueError(f"invalid oracle score for{c['id']}")
        if "analytic_optimal_distance" in c and c["analytic_optimal_distance"] != d:
            raise ValueError(f"analytic/oracle mismatch for{c['id']}")
    return data


def run_batch(args, setting, cases, repeat, batch_index, phase, scores):
    rows, batch, failures = sparse_v4.run_batch(args, setting, cases, repeat, batch_index, phase)
    batch["case_order_sha256"] = digest_json([c["id"] for c in cases])
    by_id = {c["id"]: c for c in cases}
    for row in rows:
        case = by_id[row["case_id"]]
        label = f"{case['id']} {setting} {phase}{repeat}"
        truth = scores[case["id"]]
        if row.get("distance") != truth:
            failures.append(f"{label}: distance{row.get('distance')} != independent oracle{truth}")
        for index, item in enumerate(row.get("trace", [])):
            lo, hi = item.get("lower"), item.get("upper")
            if isinstance(lo, int) and lo >= 0 and lo > truth:
                failures.append(f"{label}: round{index} lower exceeds independent optimum")
            if isinstance(hi, int) and hi >= 0 and hi < truth:
                failures.append(f"{label}: round{index} upper below independent optimum")
        row["cluster"] = case["cluster"]
        row["oracle_distance"] = truth
        row["matched_pair_id"] = case.get("matched_pair_id", "")
    return rows, batch, failures


def scopes_for(cases):
    # No all380 endpoint: controlled and random subsets have separate interpretations.
    scopes = []
    for suite in SUITES:
        chosen = [c for c in cases if c["suite"] == suite]
        scopes.append((suite, chosen))
        for n in sorted({c["source_length"] for c in chosen}):
            scopes.append((f"{suite}/n{n}", [c for c in chosen if c["source_length"] == n]))
        if suite != SUITES[0]:
            for key in ("primitive_period", "delta_percent", "regime", "direction"):
                for value in sorted({c["generation"][key] for c in chosen}):
                    scopes.append((f"{suite}/{key}={value}", [c for c in chosen if c["generation"][key] == value]))
            for n in (1024, 4096):
                for regime in ("perfect", "imperfect"):
                    scopes.append((f"{suite}/n{n}/{regime}", [c for c in chosen
                                   if c["source_length"] == n and c["generation"]["regime"] == regime]))
    return scopes


def quantile(sorted_values, p):
    index = (len(sorted_values) - 1) * p
    left = int(index)
    fraction = index - left
    return sorted_values[left] * (1 - fraction) + sorted_values[min(left + 1, len(sorted_values) - 1)] * fraction


def bootstrap_ratio(selected, lookup, baseline):
    clusters = {}
    for c in selected:
        pair = clusters.setdefault(c["cluster"], [0., 0.])
        pair[0] += lookup[(c["id"], baseline)]["metrics"]["total_ms"]["median"]
        pair[1] += lookup[(c["id"], "sparse")]["metrics"]["total_ms"]["median"]
    values = [clusters[key] for key in sorted(clusters)]
    rng = random.Random(BOOTSTRAP_SEED)
    ratios = []
    for _ in range(BOOTSTRAP_REPLICATES):
        sample = [values[rng.randrange(len(values))] for _ in values]
        ratios.append(sum(v[1] for v in sample) / sum(v[0] for v in sample))
    ratios.sort()
    return {"cluster_count": len(values), "resamples": BOOTSTRAP_REPLICATES, "seed": BOOTSTRAP_SEED,
            "percentile95_ratio_interval": [quantile(ratios, .025), quantile(ratios, .975)],
            "cluster_sizes": dict(sorted(Counter(c["cluster"] for c in selected).items())),
            "caveat": "paired synthetic-input cluster bootstrap, unadjusted subgroup intervals; not hardware variability"}


def summarize(args, cases, records, batches, failures, elapsed_ms, frozen, oracle):
    timed = [row for row in records if row["phase"] == "timed"]
    metrics = previous.numeric_keys(timed)
    by_id = {c["id"]: c for c in cases}
    groups = {}
    for row in timed:
        groups.setdefault((row["case_id"], row["bound"]), []).append(row)
    lookup, per_case = {}, []
    for (case_id, bound), rows in sorted(groups.items()):
        case = by_id[case_id]
        item = {"case_id": case_id, "bound": bound, "mode": "linear", "suite": case["suite"],
                "cluster": case["cluster"], "source_length": case["source_length"], "repetitions": len(rows), "metrics": {},
                "fallback_reasons": {key: sorted({r.get(key, "") for r in rows} - {""})
                                     for key in ("bound_fallback", "ordered_fallback", "sparse_fallback")}}
        for metric in metrics:
            values = [row[metric] for row in rows if isinstance(row.get(metric), (int, float))
                      and not isinstance(row.get(metric), bool) and math.isfinite(row[metric])]
            if values:
                item["metrics"][metric] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
        if len(rows) != args.repeats:
            failures.append(f"{case_id}/{bound}: expected7 repeats, got{len(rows)}")
        lookup[(case_id, bound)] = item
        per_case.append(item)
    expected = len(cases) * len(SETTINGS) * args.repeats
    if len(timed) != expected or len(records) != len(cases) * len(SETTINGS) * (args.repeats + 1):
        failures.append("missing timed or warmup outputs")
    if any(len({r.get("distance") for r in records if r["case_id"] == c["id"]}) != 1 for c in cases):
        failures.append("cross-setting/repetition distance disagreement")
    aggregates, comparisons, diagnostics = [], [], []
    complete = all((c["id"], b) in lookup for c in cases for b in BOUNDS)
    if complete:
        for scope, selected in scopes_for(cases):
            sums = {}
            for bound in BOUNDS:
                items = [lookup[(c["id"], bound)] for c in selected]
                a = {"scope": scope, "mode": "linear", "bound": bound, "case_count": len(items), "metrics": {}}
                for metric in metrics:
                    values = [i["metrics"][metric]["median"] for i in items if metric in i["metrics"]]
                    if values:
                        a["metrics"][metric] = {"sum_of_case_medians": sum(values),
                                                   "median_of_case_medians": statistics.median(values),
                                                   "max_of_case_medians": max(values)}
                sums[bound] = a
                aggregates.append(a)
            for baseline in BOUNDS[:-1]:
                comparison = {"scope": scope, "baseline": f"linear/{baseline}", "candidate": "linear/sparse",
                              "case_count": len(selected), "metrics": {}, "case_change_counts": {}}
                for metric in metrics:
                    left = sums[baseline]["metrics"].get(metric, {}).get("sum_of_case_medians")
                    right = sums["sparse"]["metrics"].get(metric, {}).get("sum_of_case_medians")
                    if left is not None and right is not None:
                        comparison["metrics"][metric] = {"baseline_sum_median": left, "candidate_sum_median": right,
                            "candidate_over_baseline": right / left if left else None,
                            "relative_change": (right - left) / left if left else None}
                for metric in ("total_ms", "dp_cells", "rounds", "bound_bytes"):
                    counts = {"lower": 0, "equal": 0, "higher": 0}
                    for c in selected:
                        a = lookup[(c["id"], baseline)]["metrics"][metric]["median"]
                        b = lookup[(c["id"], "sparse")]["metrics"][metric]["median"]
                        counts["lower" if b < a else "higher" if b > a else "equal"] += 1
                    comparison["case_change_counts"][metric] = counts
                comparison["bootstrap"] = bootstrap_ratio(selected, lookup, baseline)
                comparison["role"] = "co_primary" if scope == SUITES[0] and baseline in ("length", "qgram") else "secondary_descriptive"
                comparisons.append(comparison)
            diag = {"scope": scope, "case_count": len(selected), "sparse_fallback_reasons": {},
                    "same_dp_cells_as_dense": 0, "more_dp_cells_than_dense": 0, "fewer_dp_cells_than_dense": 0,
                    "sparse_fallback_cases": 0, "dense_fallback_cases": 0}
            for c in selected:
                s, d = lookup[(c["id"], "sparse")], lookup[(c["id"], "ordered")]
                sv, dv = s["metrics"]["dp_cells"]["median"], d["metrics"]["dp_cells"]["median"]
                diag["same_dp_cells_as_dense" if sv == dv else "more_dp_cells_than_dense" if sv > dv else "fewer_dp_cells_than_dense"] += 1
                reasons = s["fallback_reasons"]["sparse_fallback"]
                diag["sparse_fallback_cases"] += bool(reasons)
                diag["dense_fallback_cases"] += bool(d["fallback_reasons"]["ordered_fallback"])
                for reason in reasons:
                    diag["sparse_fallback_reasons"][reason] = diag["sparse_fallback_reasons"].get(reason, 0) + 1
            diagnostics.append(diag)
    primary = [c for c in comparisons if c["role"] == "co_primary"]
    decisions = [{"baseline": c["baseline"],
                  "candidate_over_baseline": c["metrics"]["total_ms"]["candidate_over_baseline"],
                  "percentile95_ratio_interval": c["bootstrap"]["percentile95_ratio_interval"],
                  "passes_10percent_margin": c["metrics"]["total_ms"]["candidate_over_baseline"] <= .9,
                  "passes_interval_upper": c["bootstrap"]["percentile95_ratio_interval"][1] < 1}
                 for c in primary]
    return {"schema_version": 5, "status": "fail" if failures else "pass", "frozen_configuration": frozen,
            "case_count": len(cases), "timed_result_count": len(timed), "result_count_including_warmup": len(records),
            "benchmark_wall_ms": elapsed_ms, "engine_binary_sha256": frozen["engine_binary_sha256"],
            "oracle_file_sha256": previous.sha256(args.oracle),
            "oracle_provenance": {k: v for k, v in oracle.items() if k != "scores"},
            "validation": {"failure_count": len(failures), "failures": failures,
                           "independent_oracle_case_count": len(oracle["scores"]), "analytic_reference_case_count": 160,
                           "evidence": "Every output score, certificate and finite round bound checked against independently computed full-DP optimum; every CIGAR replayed"},
            "primary_gate": {"scope": SUITES[0], "candidate": "linear/sparse", "comparisons": decisions,
                             "passed": not failures and len(decisions) == 2 and all(d["passes_10percent_margin"] and d["passes_interval_upper"] for d in decisions)},
            "timing_policy": {"primary": "sum of per-case median total_ms over7 timed repeats; all internal overhead included",
                              "overheads": "prediction, bound preprocessing, sparse proof, band DP, certificate and traceback",
                              "batch_wall": "startup and I/O separately included in batch wall_ms, not engine total_ms",
                              "selection": "one fixed candidate; no tuning or best-predictor selection",
                              "limits": "synthetic global unit-cost alignments, one runtime, no external-aligner comparison",
                              "work_proxy": "DP cells and sparse expansions are distinct units; use runtime as primary",
                              "subgroups": "descriptive unadjusted; no pooled endpoint rescues the primary gate",
                              "fallbacks": "per-case counts, not multiplied by repetitions"},
            "metric_names": metrics, "batches": batches, "per_case": per_case, "aggregates": aggregates,
            "paired_comparisons": comparisons, "sparse_diagnostics": diagnostics}


def main(argv=None):
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=root / "build" / "certalign")
    parser.add_argument("--output-dir", type=Path, default=root / "results" / "repeat_benchmark")
    parser.add_argument("--oracle", type=Path)
    parser.add_argument("--timeout", type=float, default=600.)
    parser.add_argument("--generate-only", action="store_true")
    parser.add_argument("--allow-rebuilt-binary", action="store_true",
                        help="allow a recompiled binary only in a separate output directory; freeze its actual hash")
    args = parser.parse_args(argv)
    args.binary, args.output_dir = args.binary.resolve(), args.output_dir.resolve()
    args.oracle = args.oracle.resolve() if args.oracle else args.output_dir / "oracle_scores.json"
    args.width, args.block, args.seed_k, args.repeats = 8, 64, 11, 7
    args.order_block, args.order_max_cells, args.sparse_max_work, args.bound_memory = 8, 4_000_000, 4_000_000, 64 * 1024 * 1024
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    if not args.binary.is_file():
        parser.error("engine binary does not exist")
    args.binary_sha256 = previous.sha256(args.binary)
    if args.allow_rebuilt_binary and args.output_dir == (root / "results" / "repeat_benchmark").resolve():
        parser.error("portable rebuilds require a separate --output-dir")
    if not args.allow_rebuilt_binary and args.binary_sha256 != EXPECTED_BINARY_SHA256:
        parser.error("engine binary does not match the frozen v4 binary")
    old_config = args.output_dir / "frozen_configuration.json"
    if old_config.exists() and json.loads(old_config.read_text()).get("engine_binary_sha256") != args.binary_sha256:
        parser.error("output directory already freezes a different engine binary; use a fresh output directory")
    if previous.sha256(root / "include" / "predictor.hpp") != previous.PREDICTOR_SHA256:
        parser.error("predictor source differs from fixed v1 source")
    cases = make_cases()
    frozen = freeze_inputs(args, cases, root)
    print(json.dumps({"phase": "inputs_frozen", "case_count": len(cases),
                      "configuration_sha256": frozen["configuration_sha256"],
                      "pairs_tsv_sha256": frozen["pairs_tsv_sha256"], "metadata_sha256": frozen["metadata_sha256"]}, indent=2), flush=True)
    if args.generate_only:
        return 0
    oracle = load_oracle(args.oracle, frozen, cases)
    if (args.output_dir / "benchmark_started.json").exists():
        parser.error("refusing to overwrite a started experiment; preserve the original results")
    sparse_v4.write_frozen(args.output_dir / "benchmark_started.json", json.dumps({
        "configuration_sha256": frozen["configuration_sha256"],
        "oracle_sha256": previous.sha256(args.oracle), "engine_binary_sha256": args.binary_sha256}, indent=2) + "\n")
    records, batches, failures = [], [], []
    started = time.perf_counter()
    orders = [("warmup", -1, list(SETTINGS))]
    for repeat in range(args.repeats):
        offset = repeat % len(SETTINGS)
        order = list(SETTINGS[offset:] + SETTINGS[:offset])
        if repeat % 2:
            order.reverse()
        orders.append(("timed", repeat, order))
    for phase, repeat, order in orders:
        shuffled = list(cases)
        random.Random(CASE_ORDER_SEED + repeat + 1).shuffle(shuffled)
        for setting in order:
            print(f"{phase}{repeat + 1}: {'/'.join(setting)},380 cases", file=sys.stderr, flush=True)
            rows, batch, errors = run_batch(args, setting, shuffled, repeat, len(batches), phase, oracle["scores"])
            records.extend(rows)
            batches.append(batch)
            failures.extend(errors)
            if failures:
                break
        if failures:
            break
    elapsed_ms = (time.perf_counter() - started) * 1000
    if previous.sha256(args.binary) != args.binary_sha256 or source_hashes(root) != frozen["source_sha256"]:
        failures.append("engine binary or source files changed during experiment")
    summary = summarize(args, cases, records, batches, failures, elapsed_ms, frozen, oracle)
    previous.save_results(args.output_dir, records, summary)
    print(json.dumps({"status": summary["status"], "case_count": len(cases),
                      "timed_result_count": summary["timed_result_count"], "total_outputs": len(records),
                      "failure_count": len(failures), "primary_gate": summary["primary_gate"],
                      "wall_seconds": elapsed_ms / 1000}, indent=2))
    return int(bool(failures))


if __name__ == "__main__":
    raise SystemExit(main())
