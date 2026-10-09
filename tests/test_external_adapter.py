#!/usr/bin/env python3
"""Independent score/path tests for native baseline APIs and the frozen wrapper."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import random
import subprocess
import tempfile

from test_correctness import oracle_distance, replay_cigar

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("edlib", "wfa2-high", "wfa2-ultralow", "prototype-length",
           "prototype-qgram", "prototype-sparse")


def make_cases():
    strings = ["".join(chars) for n in range(5)
               for chars in itertools.product("AC", repeat=n)]
    cases = [(f"exhaustive-{i}", x, y)
             for i, (x, y) in enumerate(itertools.product(strings, repeat=2))]
    rng = random.Random(720931)
    alphabets = ["ACGT", "ACGTNacgtn", "\x00\x01\x02\x7f\x80\xfe\xff"]
    for index in range(96):
        alphabet = alphabets[index % len(alphabets)]
        x = "".join(rng.choice(alphabet) for _ in range(rng.randrange(161)))
        y = "".join(rng.choice(alphabet) for _ in range(rng.randrange(161)))
        cases.append((f"random-{index}", x, y))
    for length in (63, 64, 65, 127, 128, 129, 255, 256, 257):
        x = "".join(rng.choice("ACGT") for _ in range(length))
        cases.extend([(f"insert-{length}", x, x[:length // 2] + "C" * 67 + x[length // 2:]),
                      (f"delete-{length}", x, x[:length // 3] + x[2 * length // 3:])])
    # Nontrivial ultralow/BiWFA work, beyond the small-alignment base cases.
    cases.extend([(f"divergent-{n}", "A" * n, "C" * n) for n in (513, 1025)])
    for index in range(4):
        x = "".join(rng.choice("ACGT") for _ in range(1025))
        y = "".join(rng.choice("ACGT") if rng.random() < .45 else c for c in x)
        y = y[:300] + "".join(rng.choice("ACGT") for _ in range(193)) + y[300:]
        cases.append((f"long-noisy-{index}", x, y))
    rng.shuffle(cases)  # Reuse after long/short/empty inputs in mixed order.
    return cases


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build/external_adapter")
    parser.add_argument("--prototype", type=Path, default=ROOT / "build/certalign")
    parser.add_argument("--output", type=Path, default=ROOT / "results/external_adapter_validation.json")
    args = parser.parse_args()
    binary = args.binary.resolve()
    prototype = args.prototype.resolve()
    cases = make_cases()
    by_id = {case[0]: case for case in cases}
    oracle = {identifier: oracle_distance(x, y) for identifier, x, y in cases}
    payload = "".join(f"{identifier}\t{x}\t{y}\n" for identifier, x, y in cases).encode("latin1")
    outputs = bounds = paths = cli_checks = 0
    native = {}
    with tempfile.TemporaryDirectory(prefix="external-adapter-test-") as directory:
        pairs = Path(directory) / "pairs.tsv"
        pairs.write_bytes(payload)
        for method in METHODS:
            completed = subprocess.run([str(binary), "--method", method, "--input", str(pairs),
                                        "--repeats", "3", "--warmups", "1"],
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            rows = [json.loads(line) for line in completed.stdout.splitlines()]
            meta, records = rows[0], rows[1:]
            assert meta["type"] == "batch" and meta["pairs"] == len(cases)
            assert len(records) == len(cases) * 3
            assert meta["batch_setup_ms"] >= 0 and meta["batch_teardown_ms"] >= 0
            assert abs(sum(row["total_ms"] for row in records) - meta["alignment_sum_ms"]) < 1e-6
            native[method] = {}
            for index, row in enumerate(records):
                identifier, x, y = cases[index % len(cases)]
                assert row["id"] == identifier and row["repeat"] == index // len(cases)
                assert row["distance"] == oracle[identifier], (method, identifier, row, oracle[identifier])
                assert row["certified_exact"] and not row["error"]
                assert replay_cigar(x, y, row["cigar"]) == row["distance"]
                paths += 1
                assert row["lower_bound"] <= oracle[identifier]
                bounds += 1
                for trace in row["trace"]:
                    if trace["lower"] >= 0:
                        assert trace["lower"] <= oracle[identifier]
                        bounds += 1
                    if trace["upper"] >= 0:
                        assert oracle[identifier] <= trace["upper"]
                        bounds += 1
                if identifier in native[method]:
                    previous = native[method][identifier]
                    assert all(row[key] == previous[key] for key in
                               ("distance", "cigar", "trace", "rounds", "dp_cells"))
                native[method][identifier] = row
                outputs += 1
        for bound in ("length", "qgram", "sparse"):
            completed = subprocess.run([str(prototype), "--mode", "linear", "--bound", bound],
                                       input=payload, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            rows = [json.loads(line) for line in completed.stdout.splitlines()]
            assert len(rows) == len(cases)
            for row in rows:
                wrapped = native["prototype-" + bound][row["id"]]
                for key in ("distance", "lower_bound", "cigar", "rounds", "dp_cells",
                            "sparse_expansions", "bound_bytes", "sparse_fallback", "certified_exact"):
                    assert row[key] == wrapped[key], (bound, row["id"], key)
                trace = [{k: (-1 if v is None else v) for k, v in round_.items()}
                         for round_ in row["trace"]]
                assert trace == wrapped["trace"]
                cli_checks += 1
    summary = {
        "passed": True, "cases": len(cases), "exhaustive_binary_pairs": 961,
        "methods": list(METHODS), "timed_repeats_in_validation": 3, "warmups": 1,
        "alignment_outputs_checked": outputs, "cigar_checks": paths,
        "finite_bound_checks": bounds, "unchanged_cli_equivalence_checks": cli_checks,
        "oracle": "independent unrestricted scalar Levenshtein DP",
        "test_pair_bytes_sha256": hashlib.sha256(payload).hexdigest(),
        "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
