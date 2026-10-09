#!/usr/bin/env python3
"""Frozen external-alignment comparison: generation, runtime freeze, run, analysis.

All workloads are semisynthetic. Nothing models a sequencer or actual read data.
Run --generate-only, independently audit/freeze, compute oracle, --freeze-runtime,
then --run. The candidate and parameters are fixed; controls are secondary.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import math
import os
import platform
import random
import re
import signal
import statistics
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/external_benchmark"
ALPHABET = "ACGT"
SEED = 2026101001
BOOTSTRAP_SEED = 2026101002
METHODS = ("edlib", "wfa2-high", "wfa2-ultralow", "prototype-length", "prototype-qgram", "prototype-sparse")
SCOPES = ("genomic", "random_control")
REPEATS = 7
SEQUENCE_SHA256 = "e999dceccb76562d40e49a10d9e6d4cf1c05d184f4f80b92bb32d4a2420c11d0"
# A fresh wrapper has exactly one child, so RUSAGE_CHILDREN is that adapter's
# peak RSS, not the maximum of earlier benchmark invocations. Linux reports KiB.
RSS_WRAPPER = ("import resource,subprocess,sys; "
               "r=subprocess.run(sys.argv[2:]); "
               "open(sys.argv[1],'w').write(str(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)); "
               "sys.exit(r.returncode)")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def hash_text(value):
    return hashlib.sha256(value.encode()).hexdigest()


def json_write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def seed_for(label):
    return int.from_bytes(hashlib.sha256(f"{SEED}/{label}".encode()).digest()[:8], "big")


def dna(rng, length):
    return "".join(rng.choice(ALPHABET) for _ in range(length))


def source_hashes():
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(
        [*ROOT.glob("include/*.hpp"), *ROOT.glob("src/*.cpp")])}


def noise_events(length, label, regime):
    rng = random.Random(seed_for(label + "/noise/" + regime))
    sub_probability, indel_probability = (.004, 0.) if regime == "low" else (.02, .002)
    deleted, inserted, substitutions = [], [], []
    for pos in range(length):
        event = rng.random() < indel_probability
        if event:
            if rng.randrange(2) == 0:
                deleted.append(pos)
                continue  # A removed base cannot also be substituted.
            inserted.append([pos, rng.choice(ALPHABET)])
        if rng.random() < sub_probability:
            substitutions.append([pos, rng.randrange(1, 4)])
    return {"substitution_probability_given_retained_base": sub_probability,
            "single_base_indel_event_probability_per_pre_noise_base": indel_probability,
            "small_deleted_positions": deleted, "small_insertions_before_position": inserted,
            "substitution_offsets": substitutions,
            "coordinate_system": "zero-based positions in target after large event, before noise",
            "operation_order": "draw small indel per pre-noise base; D removes base and skips substitution; I inserts a random base immediately before retained base; retained original base then receives Bernoulli substitution with nonzero cyclic ACGT offset; newly inserted bases receive no further noise"}


def apply_noise(sequence, events):
    deleted = set(events["small_deleted_positions"])
    inserted = dict(events["small_insertions_before_position"])
    substitutions = dict(events["substitution_offsets"])
    assert not deleted.intersection(substitutions) and not deleted.intersection(inserted)
    out = []
    for pos, base in enumerate(sequence):
        if pos in deleted:
            continue
        if pos in inserted:
            out.append(inserted[pos])
        out.append(ALPHABET[(ALPHABET.index(base) + substitutions.get(pos, 0)) % 4])
    return "".join(out)


def generate(output, reference):
    output.mkdir(parents=True, exist_ok=True)
    if (output / "raw.jsonl").exists():
        raise RuntimeError("refusing to regenerate a benchmark that already has timing results")
    fasta = gzip.decompress(reference.read_bytes()).decode()
    lines = fasta.splitlines()
    assert lines[0].startswith(">NC_000913.3 ") and sum(line.startswith(">") for line in lines) == 1
    genome = "".join(lines[1:])
    assert hash_text(genome) == SEQUENCE_SHA256 and set(genome) == set(ALPHABET)
    rng = random.Random(SEED)
    starts = []
    while len(starts) < 8:
        start = rng.randrange(len(genome) - 8192 + 1)
        if all(abs(start - previous) >= 8192 for previous in starts):
            starts.append(start)
    starts.sort()
    cases, records = [], []
    for locus, start in enumerate(starts):
        label = f"locus{locus:02d}"
        genomic_parent = genome[start:start + 8192]
        random_parent_seed = seed_for(label + "/random-parent")
        random_parent = dna(random.Random(random_parent_seed), 8192)
        for n in (1024, 4096, 8192):
            offset = (8192 - n) // 2
            parents = {"genomic": genomic_parent, "random_control": random_parent}
            for percent in (2, 8):
                gap = max(1, n * percent // 100)
                event_label = f"{label}/n{n}/p{percent}"
                event_rng = random.Random(seed_for(event_label + "/large-event"))
                position = event_rng.randrange(n // 3, n // 2 + 1)
                payload = dna(event_rng, gap)
                for direction in ("insertion", "deletion"):
                    for regime in ("low", "mixed"):
                        pair_id = f"{label}-n{n}-p{percent}-{direction}-{regime}"
                        events = noise_events(n + (gap if direction == "insertion" else -gap),
                                              event_label + "/" + direction, regime)
                        for scope in SCOPES:
                            x = parents[scope][offset:offset + n]
                            pre = x[:position] + payload + x[position:] if direction == "insertion" else x[:position] + x[position + gap:]
                            y = apply_noise(pre, events)
                            case_id = pair_id + "-" + scope
                            details = {"locus_parent_start_0based": start,
                                       "source_reference_interval_0based_half_open": [start + offset, start + offset + n] if scope == "genomic" else None,
                                       "parent_window_length": 8192, "nested_window_offset": offset,
                                       "random_parent_seed": random_parent_seed if scope == "random_control" else None,
                                       "large_event_position": position, "large_event_length": gap,
                                       "large_insertion_payload": payload if direction == "insertion" else "",
                                       "large_deleted_source": x[position:position + gap] if direction == "deletion" else "",
                                       "pre_noise_target_sha256": hash_text(pre), **events}
                            case = {"id": case_id, "scope": scope, "cluster": label,
                                    "matched_pair_id": pair_id, "source_length": n, "gap_percent": percent,
                                    "direction": direction, "noise_regime": regime,
                                    "x_length": len(x), "y_length": len(y), "generation": details,
                                    "input_sha256": hash_text(x + "\n" + y)}
                            cases.append(case)
                            records.append((case_id, x, y))
    assert len(cases) == 384 and len({c["id"] for c in cases}) == 384
    assert Counter(c["scope"] for c in cases) == {scope: 192 for scope in SCOPES}
    pairs = "".join(f"{case_id}\t{x}\t{y}\n" for case_id, x, y in records)
    (output / "pairs.tsv").write_text(pairs)
    config = {"master_seed": SEED, "bootstrap_seed": BOOTSTRAP_SEED,
              "locus_parent_starts_0based": starts, "locus_sampling": "uniform integer start with rejection of overlap; sorted after selection; no sequence/outcome filtering",
              "lengths": [1024, 4096, 8192], "gap_percents_floor": [2, 8],
              "large_event_position_range_inclusive": "floor(n/3)..floor(n/2)",
              "noise_regimes": ["low", "mixed"], "cases": 384, "primary_scope": "genomic",
              "control_scope": "random_control", "matched_masks": "Within each locus/length/gap/direction/regime, genomic and random use identical large-event position, insertion payload, small-indel positions/directions/bases and nonzero substitution offsets. Deletion removes actual source bases.",
              "noise_note": "low: .004 substitutions; mixed: .02 substitutions plus .002 single-base I/D event probability with equiprobable direction. Substitution probability conditional on retained pre-noise base.",
              "evidence_status": "semisynthetic reference-derived sequences with planted errors; not actual reads, not platform-calibrated",
              "methods": METHODS, "candidate": "prototype-sparse", "timed_repeats": REPEATS,
              "warmup_repeats": 1, "bootstrap_replicates": 9999,
              "primary_success_rule": "All outputs correct, no candidate failure, candidate aggregate runtime <=0.9 of each of edlib, wfa2-high and wfa2-ultralow; each paired locus-cluster bootstrap 95% upper runtime ratio <1. No per-case baseline selection.",
              "prototype_parameters": {"mode": "linear", "width": 8, "max_cells": 50000000,
                  "epsilon_num": 0, "epsilon_den": 1, "bound_memory": 67108864, "order_block": 8, "sparse_max_work": 4000000},
              "reference_accession": "NC_000913.3", "reference_sequence_sha256": SEQUENCE_SHA256,
              "reference_compressed_sha256": sha(reference)}
    json_write(output / "config.json", config)
    json_write(output / "cases.json", {"cases": cases})
    frozen = {"stage": "inputs frozen before candidate alignment or performance timing", "source_sha256": source_hashes(),
              "files": {name: sha(output / name) for name in ("pairs.tsv", "cases.json", "config.json")},
              "generator_sha256": sha(__file__)}
    json_write(output / "input_freeze.json", frozen)
    print(json.dumps({"status": "generated", "count": len(cases), "input_freeze": frozen, "starts": starts}, indent=2))


def verify_inputs(output):
    frozen = json.loads((output / "input_freeze.json").read_text())
    for name, digest in frozen["files"].items():
        assert sha(output / name) == digest, f"input changed: {name}"
    assert source_hashes() == frozen["source_sha256"], "v5 engine changed"
    cases = json.loads((output / "cases.json").read_text())["cases"]
    records = {}
    for line in (output / "pairs.tsv").read_text().splitlines():
        case_id, x, y = line.split("\t")
        assert case_id not in records
        records[case_id] = (x, y)
    assert len(cases) == len(records) == 384
    for case in cases:
        x, y = records[case["id"]]
        assert hash_text(x + "\n" + y) == case["input_sha256"]
    return cases, records


def freeze_runtime(output, adapter):
    verify_inputs(output)
    if (output / "raw.jsonl").exists():
        raise RuntimeError("cannot freeze after timing begins")
    paths = {"harness": Path(__file__), "adapter": adapter,
             "build_manifest": ROOT / "build/external_build.json",
             "adapter_source": ROOT / "benchmarks/external_adapter.cpp",
             "build_script": ROOT / "benchmarks/build_external.py",
             "plan": ROOT / "EXTERNAL_BENCHMARK_PLAN.md",
             "oracle_scores": output / "oracle_scores.json"}
    manifest = {"files": {name: {"path": str(path.relative_to(ROOT)), "sha256": sha(path)} for name, path in paths.items()},
                "input_freeze_sha256": sha(output / "input_freeze.json"),
                "per_batch_timeout_seconds": 300,
                "source_sha256": source_hashes(), "platform": platform.platform(),
                "python": sys.version, "timing_started": False}
    json_write(output / "runtime_freeze.json", manifest)
    # Preserve the build environment and dependency revisions as a committed result.
    json_write(output / "build_manifest.json", json.loads(paths["build_manifest"].read_text()))
    print(json.dumps({"status": "runtime frozen", "files": manifest["files"]}, indent=2))


def verify_runtime(output):
    verify_inputs(output)
    manifest = json.loads((output / "runtime_freeze.json").read_text())
    assert sha(output / "input_freeze.json") == manifest["input_freeze_sha256"]
    for name, spec in manifest["files"].items():
        assert sha(ROOT / spec["path"]) == spec["sha256"], f"runtime freeze changed: {name}"
    assert source_hashes() == manifest["source_sha256"]


CIGAR = re.compile(r"([1-9][0-9]*)([=XID])")


def validate(row, case, pair, optimum):
    failures = []
    if row.get("error"):
        failures.append("alignment error: " + row["error"])
    if row.get("distance") != optimum:
        failures.append(f"distance {row.get('distance')} != oracle {optimum}")
    if row.get("certified_exact") is not True:
        failures.append("not certified exact")
    x, y = pair
    if row.get("n") != len(x) or row.get("m") != len(y):
        failures.append("input length mismatch")
    lower = row.get("lower_bound", -1)
    if lower >= 0 and lower > optimum:
        failures.append("lower bound exceeds oracle")
    for index, round_trace in enumerate(row.get("trace", [])):
        round_lower = round_trace.get("lower", -1)
        round_upper = round_trace.get("upper", -1)
        if round_lower >= 0 and round_lower > optimum:
            failures.append(f"round {index} lower bound exceeds oracle")
        if round_upper >= 0 and round_upper < optimum:
            failures.append(f"round {index} upper bound below oracle")
    cigar = row.get("cigar", "")
    i = j = distance = end = 0
    try:
        for match in CIGAR.finditer(cigar):
            assert match.start() == end
            count, op = int(match[1]), match[2]
            end = match.end()
            if op in "=X":
                assert i + count <= len(x) and j + count <= len(y)
                assert all((a == b) == (op == "=") for a, b in zip(x[i:i + count], y[j:j + count]))
                i += count
                j += count
                distance += count * (op == "X")
            elif op == "D":
                i += count
                distance += count
            else:
                j += count
                distance += count
        assert end == len(cigar) and i == len(x) and j == len(y) and distance == optimum
    except (AssertionError, ValueError):
        failures.append("invalid or nonoptimal CIGAR")
    if not isinstance(row.get("total_ms"), (float, int)) or not math.isfinite(row["total_ms"]) or row["total_ms"] < 0:
        failures.append("invalid runtime")
    return failures


def run(output, adapter):
    verify_runtime(output)
    cases, records = verify_inputs(output)
    oracle = json.loads((output / "oracle_scores.json").read_text())
    assert oracle["pairs_sha256"] == sha(output / "pairs.tsv") and oracle["metadata_sha256"] == sha(output / "cases.json")
    assert oracle["case_count"] == 384 and set(oracle["scores"]) == set(records)
    if (output / "raw.jsonl").exists():
        raise RuntimeError("timing results already exist; refusing unplanned rerun")
    by_id = {c["id"]: c for c in cases}
    schedule = []
    for repeat in range(-1, REPEATS):
        methods = list(METHODS)
        shift = (repeat + 1) % len(methods)
        methods = methods[shift:] + methods[:shift]
        if repeat % 2:
            methods.reverse()
        scopes = SCOPES if repeat % 2 == 0 else tuple(reversed(SCOPES))
        for scope in scopes:
            ids = [c["id"] for c in cases if c["scope"] == scope]
            random.Random(seed_for(f"case-order/{repeat}/{scope}")).shuffle(ids)
            path = output / f"order-{repeat + 1}-{scope}.tsv"
            path.write_text("".join(f"{i}\t{records[i][0]}\t{records[i][1]}\n" for i in ids))
            for method in methods:
                schedule.append({"repeat": repeat, "scope": scope, "method": method,
                                 "input": path.name, "input_sha256": sha(path), "ids": ids})
    json_write(output / "schedule.json", schedule)
    failures, batch_count, row_count, finite_bounds, finite_round_lower, finite_round_upper = [], 0, 0, 0, 0, 0
    started = time.time()
    with (output / "raw.jsonl").open("w") as raw:
        for task in schedule:
            command = [str(adapter), "--method", task["method"], "--input", str(output / task["input"]), "--repeats", "1", "--warmups", "0"]
            rss_path = output / ".last_process_rss"
            timeout_failure = None
            process = subprocess.Popen([sys.executable, "-c", RSS_WRAPPER, str(rss_path), *command],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                       start_new_session=True)
            try:
                stdout, stderr = process.communicate(timeout=300)
            except subprocess.TimeoutExpired:
                # Kill the entire group, including the RSS wrapper's adapter child.
                # A surviving timed-out child must not contaminate later batches.
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                stdout, stderr = process.communicate()
                timeout_failure = "adapter exceeded fixed 300-second per-batch timeout"
            completed = subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
            try:
                process_rss_kib = int(float(rss_path.read_text().splitlines()[-1]))
            except (OSError, ValueError, IndexError):
                process_rss_kib = None
            rss_path.unlink(missing_ok=True)
            try:
                outputs = [json.loads(line) for line in completed.stdout.splitlines()]
            except json.JSONDecodeError:
                outputs = []
            batches = [row for row in outputs if row.get("type") == "batch"]
            rows = [row for row in outputs if row.get("type") == "alignment"]
            invocation_error = timeout_failure
            if completed.returncode != 0 or len(batches) != 1 or [row.get("id") for row in rows] != task["ids"]:
                invocation_error = invocation_error or f"adapter invocation failed/incomplete: exit {completed.returncode}; batches {len(batches)}; rows {len(rows)}"
            for batch in batches:
                valid_batch = (batch.get("method") == task["method"] and batch.get("pairs") == 192
                               and batch.get("repeats") == 1 and batch.get("warmups") == 0)
                for key in ("batch_setup_ms", "batch_teardown_ms", "alignment_sum_ms", "wallcompute_ms", "warmup_alignment_ms"):
                    value = batch.get(key)
                    valid_batch = valid_batch and isinstance(value, (int, float)) and math.isfinite(value) and value >= 0
                if not valid_batch:
                    invocation_error = invocation_error or "invalid batch metadata/timing"
            context = {"outer_repeat": task["repeat"], "scope": task["scope"], "warmup": task["repeat"] < 0}
            for batch in batches:
                batch.update(context)
                batch["returncode"] = completed.returncode
                batch["peak_process_rss_kib"] = process_rss_kib
                raw.write(json.dumps(batch, sort_keys=True) + "\n")
                batch_count += 1
            seen = set()
            for row in rows:
                row.update(context)
                case_id = row.get("id")
                if case_id not in by_id or case_id in seen:
                    row["validation_errors"] = ["unexpected/duplicate case id"]
                else:
                    seen.add(case_id)
                    row["validation_errors"] = validate(row, by_id[case_id], records[case_id], oracle["scores"][case_id])
                    finite_bounds += row.get("lower_bound", -1) >= 0
                    finite_round_lower += sum(t.get("lower", -1) >= 0 for t in row.get("trace", []))
                    finite_round_upper += sum(t.get("upper", -1) >= 0 for t in row.get("trace", []))
                if row.get("method") != task["method"]:
                    row["validation_errors"].append("row method differs from invoked method")
                if row["validation_errors"]:
                    failures.append({"method": task["method"], "id": case_id, **context, "errors": row["validation_errors"]})
                raw.write(json.dumps(row, sort_keys=True) + "\n")
                row_count += 1
            if invocation_error:
                failure = {"type": "invocation_failure", "method": task["method"], **context,
                           "error": invocation_error, "stderr": completed.stderr,
                           "missing_ids": sorted(set(task["ids"]) - seen)}
                raw.write(json.dumps(failure, sort_keys=True) + "\n")
                failures.append(failure)
            raw.flush()
            print(f"{task['repeat']:2} {task['scope']:14} {task['method']:18}: {len(rows)} outputs, cumulative failures {len(failures)}", flush=True)
    verify_runtime(output)
    json_write(output / "validation.json", {"all_correct": not failures,
        "expected_alignment_outputs": 384 * len(METHODS) * (REPEATS + 1),
        "alignment_outputs": row_count, "batch_outputs": batch_count, "finite_lower_bounds_checked": finite_bounds,
        "finite_round_lower_bounds_checked": finite_round_lower, "finite_round_upper_bounds_checked": finite_round_upper,
        "oracle_scores": 384, "failures": failures,
        "elapsed_seconds_including_subprocess_io_validation": time.time() - started})
    analyze(output)


def percentile(values, probability):
    position = (len(values) - 1) * probability
    lower = math.floor(position)
    fraction = position - lower
    return values[lower] * (1 - fraction) + values[min(lower + 1, len(values) - 1)] * fraction


def analyze(output):
    cases, _ = verify_inputs(output)
    by_id = {c["id"]: c for c in cases}
    rows, batches = defaultdict(list), defaultdict(list)
    all_rows = [json.loads(line) for line in (output / "raw.jsonl").read_text().splitlines()]
    for row in all_rows:
        if row.get("warmup"):
            continue
        if row["type"] == "alignment":
            rows[(row["method"], row["id"])].append(row)
        elif row["type"] == "batch":
            batches[(row["method"], row["scope"])].append(row)
    validation = json.loads((output / "validation.json").read_text())
    medians, aggregate = {}, {}
    per_case = []
    for scope in SCOPES:
        ids = [c["id"] for c in cases if c["scope"] == scope]
        for method in METHODS:
            batch = batches[(method, scope)]
            failures = [f for f in validation["failures"] if f.get("method") == method and f.get("scope") == scope]
            complete = not failures and len(batch) == REPEATS and all(len(rows[(method, i)]) == REPEATS for i in ids)
            if not complete:
                aggregate[f"{scope}/{method}"] = {"complete_and_correct": False, "total_ms": None, "failure_count": len(failures)}
                continue
            overhead = statistics.median(b["batch_setup_ms"] + b["batch_teardown_ms"] for b in batch)
            values = {}
            numeric_keys = set.intersection(*(set(k for k, v in r.items() if isinstance(v, (int, float)) and not isinstance(v, bool)) for i in ids for r in rows[(method, i)]))
            metric_sums = Counter()
            for case_id in ids:
                samples = rows[(method, case_id)]
                metrics = {k: statistics.median(r[k] for r in samples) for k in sorted(numeric_keys)}
                metric_sums.update(metrics)
                values[case_id] = metrics["total_ms"] + overhead / len(ids)
                per_case.append({"id": case_id, "scope": scope, "method": method,
                                 "median": metrics, "amortized_setup_teardown_ms": overhead / len(ids),
                                 "adjusted_median_ms": values[case_id],
                                 "sparse_fallback_values": dict(Counter(r.get("sparse_fallback", "") for r in samples))})
            medians[(method, scope)] = values
            actual_totals = [b["alignment_sum_ms"] + b["batch_setup_ms"] + b["batch_teardown_ms"] for b in batch]
            aggregate[f"{scope}/{method}"] = {"complete_and_correct": True,
                "total_ms": sum(values.values()), "sum_per_case_median_ms": metric_sums["total_ms"],
                "median_setup_teardown_ms": overhead, "median_actual_batch_total_ms": statistics.median(actual_totals),
                "actual_batch_totals_ms": actual_totals,
                "actual_batch_wallcompute_ms": [b["wallcompute_ms"] for b in batch],
                "median_actual_batch_wallcompute_ms": statistics.median(b["wallcompute_ms"] for b in batch),
                "sum_median_metrics": dict(metric_sums),
                "peak_process_rss_kib": max((b.get("peak_process_rss_kib") or 0 for b in batch), default=0),
                "rss_scope": "whole adapter process including preloaded input, all retained results and library workspaces; not per-alignment auxiliary memory"}
    comparisons = []
    rng = random.Random(BOOTSTRAP_SEED)
    clusters = sorted({c["cluster"] for c in cases})
    bootstrap = [[rng.choice(clusters) for _ in clusters] for _ in range(9999)]
    groups = {scope: [c["id"] for c in cases if c["scope"] == scope] for scope in SCOPES}
    for scope in SCOPES:
        for field in ("source_length", "gap_percent", "direction", "noise_regime"):
            for value in sorted({c[field] for c in cases}, key=str):
                groups[f"{scope}/{field}={value}"] = [c["id"] for c in cases if c["scope"] == scope and c[field] == value]
    for group, ids in groups.items():
        scope = by_id[ids[0]]["scope"]
        candidate = medians.get(("prototype-sparse", scope))
        for baseline in METHODS[:-1]:
            values = medians.get((baseline, scope))
            if candidate is None or values is None:
                comparisons.append({"group": group, "baseline": baseline, "valid": False})
                continue
            cand_cluster, base_cluster = Counter(), Counter()
            for case_id in ids:
                cluster = by_id[case_id]["cluster"]
                cand_cluster[cluster] += candidate[case_id]
                base_cluster[cluster] += values[case_id]
            ratios = sorted(sum(cand_cluster[c] for c in sample) / sum(base_cluster[c] for c in sample) for sample in bootstrap)
            ratio = sum(candidate[i] for i in ids) / sum(values[i] for i in ids)
            comparisons.append({"group": group, "baseline": baseline, "valid": True, "cases": len(ids),
                "candidate_total_ms": sum(candidate[i] for i in ids), "baseline_total_ms": sum(values[i] for i in ids),
                "runtime_ratio": ratio, "change_percent": 100 * (ratio - 1),
                "paired_locus_bootstrap_95_percentile_ci": [percentile(ratios, .025), percentile(ratios, .975)],
                "cluster_count": len(clusters), "evidence": "primary" if group == "genomic" and baseline in METHODS[:3] else "secondary diagnostic"})
    primary = [c for c in comparisons if c["group"] == "genomic" and c["baseline"] in METHODS[:3]]
    success = validation["all_correct"] and all(c["valid"] and c["runtime_ratio"] <= .9 and c["paired_locus_bootstrap_95_percentile_ci"][1] < 1 for c in primary)
    summary = {"primary_success": success, "primary_scope": "genomic", "candidate": "prototype-sparse",
        "all_outputs_correct": validation["all_correct"], "aggregation": "sum of each case's median complete native alignment time over 7 repetitions plus median per-scope batch setup+teardown; bootstrap allocates that overhead equally across the 192 cases",
        "batch_scope": "each invocation contains exactly one 192-case scope, one repetition; WFA workspace setup/teardown shared within that scope; all I/O and CIGAR verification outside measured native scope",
        "warmup": "one discarded 192-case invocation for each method and scope before all timed repetitions; timed invocations create fresh WFA workspaces",
        "bootstrap": {"seed": BOOTSTRAP_SEED, "replicates": 9999, "clusters": 8,
            "unit": "reference locus with nested lengths and all event/noise variants retained together; matched controls share cluster", "limitation": "eight independent genomic windows; CI does not represent platform, genome or hardware generalization"},
        "aggregate": aggregate, "comparisons": comparisons,
        "raw_sha256": sha(output / "raw.jsonl"), "failures_not_credited_as_speedups": True}
    json_write(output / "per_case.json", per_case)
    json_write(output / "summary.json", summary)
    print(json.dumps({"primary_success": success, "primary_comparisons": primary}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--reference", type=Path, default=ROOT / "data/reference/NC_000913.3.fna.gz")
    parser.add_argument("--adapter", type=Path, default=ROOT / "build/external_adapter")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--generate-only", action="store_true")
    mode.add_argument("--freeze-runtime", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--analyze", action="store_true")
    args = parser.parse_args()
    if args.generate_only:
        generate(args.output, args.reference)
    elif args.freeze_runtime:
        freeze_runtime(args.output, args.adapter)
    elif args.run:
        run(args.output, args.adapter)
    else:
        analyze(args.output)


if __name__ == "__main__":
    main()
