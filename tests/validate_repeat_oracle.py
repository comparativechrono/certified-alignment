#!/usr/bin/env python3
"""Build independent exact references for the frozen fresh-repeat experiment.

The helper computes unrestricted two-row Levenshtein DP. It never imports the
alignment engine or reads its results. A separate Python oracle first checks
the compiled helper on exhaustive short strings and randomized literal bytes.
The large-input oracle pass is untimed with respect to algorithm benchmarking.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import time

from test_correctness import oracle_distance

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests" / "reference_distance.cpp"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_oracle(binary: Path, records: list[tuple[str, bytes, bytes]]) -> dict[str, int]:
    payload = b"".join(case_id.encode("ascii") + b"\t" + x + b"\t" + y + b"\n"
                       for case_id, x, y in records)
    completed = subprocess.run([str(binary)], input=payload, capture_output=True,
                               check=True)
    lines = completed.stdout.splitlines()
    assert len(lines) == len(records), (len(lines), len(records))
    result: dict[str, int] = {}
    for record, line in zip(records, lines):
        output = json.loads(line)
        case_id = record[0]
        assert output["id"] == case_id, (case_id, output)
        value = output["distance"]
        assert isinstance(value, int) and not isinstance(value, bool) and value >= 0
        assert case_id not in result, f"duplicate id {case_id}"
        result[case_id] = value
    return result


def self_test(binary: Path) -> dict[str, int]:
    words = [bytes(chars) for n in range(5)
             for chars in itertools.product(b"AC", repeat=n)]
    records = [(f"exhaustive_{i}_{j}", x, y)
               for i, x in enumerate(words) for j, y in enumerate(words)]
    exhaustive_count = len(records)
    rng = random.Random(709731)
    # Latin-1 gives a one-to-one byte-to-codepoint mapping to the Python oracle.
    # Only the TSV field and record separators are excluded from byte inputs.
    byte_alphabet = bytes(value for value in range(256) if value not in (9, 10))
    for label, alphabet in (("byte", byte_alphabet), ("dna", b"ACGT")):
        for index in range(100):
            x = bytes(rng.choice(alphabet) for _ in range(rng.randrange(201)))
            y = bytes(rng.choice(alphabet) for _ in range(rng.randrange(201)))
            records.append((f"random_{label}_{index}", x, y))
    records.extend([
        ("empty_both", b"", b""),
        ("empty_left", b"", byte_alphabet),
        ("empty_right", byte_alphabet, b""),
        ("all_bytes_identity", byte_alphabet, byte_alphabet),
        ("all_bytes_reverse", byte_alphabet, byte_alphabet[::-1]),
        ("long_identity", b"ACGT" * 250, b"ACGT" * 250),
        ("long_insertion", b"AC" * 200, b"AC" * 100 + b"T" * 71 + b"AC" * 100),
    ])
    observed = run_oracle(binary, records)
    for case_id, x, y in records:
        expected = oracle_distance(x.decode("latin1"), y.decode("latin1"))
        assert observed[case_id] == expected, (case_id, observed[case_id], expected)
    for malformed in (b"no-tabs\n", b"only\tone\n", b"too\tmany\ttabs\there\n"):
        response = subprocess.run([str(binary)], input=malformed, capture_output=True)
        assert response.returncode != 0, "malformed TSV accepted"
    return {"exhaustive_binary_length_max": 4,
            "exhaustive_pairs": exhaustive_count,
            "random_byte_pairs": 100,
            "random_dna_pairs": 100,
            "additional_edge_pairs": 7,
            "python_crosschecks": len(records),
            "malformed_input_checks": 3}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "reference_distance")
    parser.add_argument("--pairs", type=Path,
                        default=ROOT / "results" / "repeat_benchmark" / "pairs.tsv")
    parser.add_argument("--metadata", type=Path,
                        default=ROOT / "results" / "repeat_benchmark" / "cases.json")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "repeat_benchmark" / "oracle_scores.json")
    parser.add_argument("--expected-count", type=int, default=380)
    parser.add_argument("--expected-analytic", type=int, default=160)
    parser.add_argument("--self-test-only", action="store_true")
    args = parser.parse_args()
    binary = args.binary.resolve()
    before_binary_sha = sha256(binary)
    before_source_sha = sha256(SOURCE)
    self_validation = self_test(binary)
    if args.self_test_only:
        print(json.dumps({"status": "passed", "validation": self_validation}, indent=2))
        return

    metadata_bytes = args.metadata.read_bytes()
    pairs_bytes = args.pairs.read_bytes()
    metadata = json.loads(metadata_bytes)
    cases = metadata["cases"]
    records = []
    for line in pairs_bytes.split(b"\n"):
        if not line:
            continue
        fields = line.split(b"\t")
        assert len(fields) == 3, "malformed dataset TSV"
        case_id = fields[0].decode("ascii")
        records.append((case_id, fields[1], fields[2]))
    assert len(records) == args.expected_count, (len(records), args.expected_count)
    assert len(cases) == len(records)
    metadata_by_id = {case["id"]: case for case in cases}
    assert len(metadata_by_id) == len(cases), "duplicate metadata id"
    assert {record[0] for record in records} == set(metadata_by_id)
    for case_id, x, y in records:
        case = metadata_by_id[case_id]
        assert len(x) == case["x_length"] and len(y) == case["y_length"], case_id
        assert hashlib.sha256(x + b"\n" + y).hexdigest() == case["input_sha256"], case_id
        assert set(x + y) <= set(b"ACGT"), f"non-DNA byte in {case_id}"
    started = time.perf_counter()
    scores = run_oracle(binary, records)
    oracle_elapsed = time.perf_counter() - started
    analytic = 0
    for case in cases:
        expected = case.get("analytic_optimal_distance")
        if expected is not None:
            assert scores[case["id"]] == expected, (case["id"], scores[case["id"]], expected)
            analytic += 1
    assert analytic == args.expected_analytic, (analytic, args.expected_analytic)
    assert pairs_bytes == args.pairs.read_bytes(), "pairs changed during reference pass"
    assert metadata_bytes == args.metadata.read_bytes(), "metadata changed during reference pass"
    assert before_binary_sha == sha256(binary), "oracle binary changed during reference pass"
    assert before_source_sha == sha256(SOURCE), "oracle source changed during reference pass"
    result = {
        "schema_version": 1,
        "method": "independent unrestricted two-row unit-cost Levenshtein DP",
        "benchmark_timing_includes_oracle": False,
        "pairs_sha256": hashlib.sha256(pairs_bytes).hexdigest(),
        "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
        "oracle_source_sha256": before_source_sha,
        "oracle_binary_sha256": before_binary_sha,
        "case_count": len(records),
        "total_dp_cells": sum(len(x) * len(y) for _, x, y in records),
        "dp_cells_definition": "interior recurrence cells n*m, excluding boundary initialization",
        "oracle_elapsed_seconds": oracle_elapsed,
        "scores": scores,
        "validation": {**self_validation, "analytic_distance_checks": analytic,
                       "dataset_input_hash_checks": len(records)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({key: value for key, value in result.items() if key != "scores"}, indent=2))


if __name__ == "__main__":
    main()
