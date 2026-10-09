#!/usr/bin/env python3
"""Independent end-to-end tests for the certified alignment prototype.

Uses an ordinary full-matrix edit-distance oracle, never the engine's band
construction, stopping rule, anchor selection, or lower-bound calculation.
All CLI requests are batched through stdin. Run after building certalign.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path
import random
import re
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
MODES = ("full", "linear", "cost", "anchor-linear", "anchor-cost")
TRACEBACK_MODES = set(MODES)
CIGAR_TOKEN = re.compile(r"([1-9][0-9]*)([=XID])")


def oracle_distance(x: str, y: str) -> int:
    """Unrestricted Levenshtein DP with unit insertion/deletion/substitution."""
    previous = list(range(len(y) + 1))
    for i, a in enumerate(x, 1):
        current = [i]
        for j, b in enumerate(y, 1):
            current.append(min(previous[j] + 1,
                               current[j - 1] + 1,
                               previous[j - 1] + (a != b)))
        previous = current
    return previous[-1]


def replay_cigar(x: str, y: str, cigar: str) -> int:
    """Check the returned witness independently, including = versus X."""
    assert isinstance(cigar, str), "CIGAR must be a string"
    if not cigar:
        assert not x and not y, "empty CIGAR for nonempty input"
        return 0
    i = j = score = position = 0
    for match in CIGAR_TOKEN.finditer(cigar):
        assert match.start() == position, f"invalid CIGAR at offset {position}"
        count, operation = int(match.group(1)), match.group(2)
        position = match.end()
        if operation in "=X":
            assert i + count <= len(x) and j + count <= len(y), "CIGAR overruns input"
            a, b = x[i:i + count], y[j:j + count]
            if operation == "=":
                assert a == b, "= operation contains unequal symbols"
            else:
                assert all(u != v for u, v in zip(a, b)), "X operation contains an equal pair"
                score += count
            i += count
            j += count
        elif operation == "D":
            i += count
            score += count
        else:
            j += count
            score += count
        assert i <= len(x) and j <= len(y), "CIGAR overruns input"
    assert position == len(cigar), "invalid trailing CIGAR"
    assert (i, j) == (len(x), len(y)), "CIGAR does not consume both inputs"
    return score


@dataclass(frozen=True)
class Case:
    name: str
    x: str
    y: str
    truth: int


def make_cases() -> tuple[list[Case], list[Case]]:
    alphabet = "AC"
    words = ["".join(t) for n in range(5)
             for t in itertools.product(alphabet, repeat=n)]
    exhaustive = [Case(f"binary_{i}_{j}", x, y, oracle_distance(x, y))
                  for i, x in enumerate(words) for j, y in enumerate(words)]
    pairs: list[tuple[str, str, str]] = [
        ("both_empty", "", ""),
        ("empty_left", "", "ACGT" * 29),
        ("empty_right", "ACGT" * 29, ""),
        ("identical_repeat", "AC" * 90, "AC" * 90),
        ("all_mismatch", "A" * 97, "T" * 101),
        ("homopolymer_imbalance", "A" * 3, "A" * 151),
        ("repeat_rotation", "ACGT" * 40, "GTAC" * 40),
        ("opposed_runs", "A" * 45 + "C" * 45, "C" * 45 + "A" * 45),
        ("compensating_indels", "A" * 36 + "G" * 70 + "T" * 36,
         "C" * 36 + "A" * 36 + "G" * 70),
        ("long_insertion", "ACGT" * 20 + "GATTACA" * 9,
         "ACGT" * 20 + "T" * 85 + "GATTACA" * 9),
        ("long_deletion", "ACGT" * 20 + "T" * 85 + "GATTACA" * 9,
         "ACGT" * 20 + "GATTACA" * 9),
        ("repeat_false_anchors", "ACGT" * 24 + "A" * 29 + "TGCA" * 19,
         "TGCA" * 19 + "A" * 29 + "ACGT" * 24),
        ("reversed_unique_chunks", "ACGTTGCAACGATCGTACCT" * 3 + "G" * 37,
         "G" * 37 + "ACGTTGCAACGATCGTACCT" * 3),
        ("single_substitution", "A" * 71 + "C" + "G" * 68,
         "A" * 71 + "T" + "G" * 68),
        ("epsilon_boundary", "ACAC", "CACA"),
    ]
    rng = random.Random(731991)
    # Longer than the predictor's anchor spacing (4 * block = 256), so these
    # exercise actual interior anchors rather than only endpoint fallbacks.
    long_unique = "".join(rng.choice("ACGT") for _ in range(820))
    pairs.extend([
        ("anchor_long_insertion", long_unique,
         long_unique[:310] + "TGCATGCA" * 13 + long_unique[310:]),
        ("anchor_long_deletion", long_unique,
         long_unique[:340] + long_unique[445:]),
        ("anchor_compensating_indels", long_unique,
         long_unique[:290] + "T" * 83 + long_unique[290:560] + long_unique[643:]),
        ("anchor_reordered_unique", long_unique,
         long_unique[:270] + long_unique[550:790] + long_unique[270:550] + long_unique[790:]),
        ("non_dna_symbols", "aN?7-.*" * 17, "aN?7-.*" * 8 + "9#[]" + "aN?7-.*" * 9),
        ("iupac_and_case", "ACGTNRYSwkmBDHV" * 11, "acgtNRYSwkmBDHV" * 11),
    ])
    for index in range(120):
        x = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 151)))
        y = list(x)
        for _ in range(rng.randrange(1, 35)):
            action = rng.randrange(3)
            if action == 0 and y:
                del y[rng.randrange(len(y))]
            elif action == 1 and y:
                y[rng.randrange(len(y))] = rng.choice("ACGT")
            else:
                at = rng.randrange(len(y) + 1)
                y[at:at] = [rng.choice("ACGT") for _ in range(rng.randrange(1, 5))]
        pairs.append((f"mutation_{index}", x, "".join(y)))
    for index in range(30):
        x = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 81)))
        y = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 81)))
        pairs.append((f"independent_{index}", x, y))
    stress = [Case(name, x, y, oracle_distance(x, y)) for name, x, y in pairs]
    return exhaustive, stress


@dataclass
class Audit:
    binary: Path
    checks: int = 0
    finite_trace_rounds: int = 0
    null_trace_rounds: int = 0
    cigar_checks: int = 0
    exact_results: int = 0
    approximate_results: int = 0
    budget_errors: int = 0
    batches: list[dict] = field(default_factory=list)
    failures: list[dict] = field(default_factory=list)

    @staticmethod
    def finite(value: object) -> bool:
        return (isinstance(value, (int, float)) and not isinstance(value, bool)
                and math.isfinite(value))

    def verify(self, result: dict, case: Case, mode: str,
               epsilon: str, allow_budget_error: bool) -> None:
        assert result.get("id") == case.name, "output id/order does not match input"
        truth = case.truth
        trace = result.get("trace", [])
        assert isinstance(trace, list), "trace is not a list"
        for index, record in enumerate(trace):
            upper, lower = record.get("upper"), record.get("lower")
            if self.finite(lower):
                assert 0 <= lower <= truth, f"round {index}: lower={lower} > true={truth}"
            if self.finite(upper):
                self.finite_trace_rounds += 1
                assert upper >= truth, f"round {index}: upper={upper} < true={truth}"
                if self.finite(lower):
                    assert lower <= upper, f"round {index}: reversed bounds"
            else:
                self.null_trace_rounds += 1
                assert upper is None, f"round {index}: nonfinite upper must be JSON null"
            if lower is not None:
                assert self.finite(lower), f"round {index}: nonfinite lower must be JSON null"
            if "cells" in record:
                assert isinstance(record["cells"], int) and record["cells"] >= 0
        error = result.get("error")
        if error:
            assert allow_budget_error, f"unexpected engine error: {error}"
            assert not result.get("certified_exact", False), "error falsely certifies exactness"
            assert not result.get("certified_approx", False), "error falsely certifies approximation"
            self.budget_errors += 1
            return
        distance = result.get("distance")
        lower = result.get("lower_bound")
        assert self.finite(distance), "successful result has no finite distance"
        assert self.finite(lower), "successful result has no finite lower bound"
        assert int(distance) == distance and int(lower) == lower, "unit-cost bounds must be integers"
        assert 0 <= lower <= truth <= distance, (
            f"final interval [{lower}, {distance}] excludes true distance {truth}")
        exact = result.get("certified_exact", False)
        if exact:
            self.exact_results += 1
            assert distance == truth, f"false exact certificate: {distance} != {truth}"
        accuracy = Fraction(epsilon)
        if accuracy == 0:
            assert exact, "zero-epsilon run terminated without exact certificate"
        else:
            assert result.get("certified_approx", False) or exact, "missing approximation certificate"
            # Compare integers against the exact decimal/rational request; a
            # floating rounding tolerance must not excuse a violated guarantee.
            assert distance * accuracy.denominator <= truth * (
                accuracy.numerator + accuracy.denominator), (
                f"distance {distance} violates (1+{epsilon}) * optimum {truth}")
            if not exact:
                self.approximate_results += 1
                assert distance * accuracy.denominator <= lower * (
                    accuracy.numerator + accuracy.denominator), (
                    "approximation flag unsupported by reported final lower bound")
        if mode in TRACEBACK_MODES:
            score = replay_cigar(case.x, case.y, result.get("cigar"))
            assert score == distance, f"CIGAR score {score} != reported distance {distance}"
            self.cigar_checks += 1
        if "dp_cells" in result:
            assert isinstance(result["dp_cells"], int) and result["dp_cells"] >= 0
            assert result["dp_cells"] >= max((t.get("cells", 0) for t in trace), default=0)

    def batch(self, name: str, cases: list[Case], mode: str, *, width: int = 0,
              epsilon: str = "0", max_cells: int = 100_000_000,
              allow_budget_error: bool = False) -> None:
        command = [str(self.binary), "--mode", mode, "--width", str(width),
                   "--block", "64", "--seed-k", "11", "--epsilon", epsilon,
                   "--max-cells", str(max_cells)]
        request = "".join(f"{c.name}\t{c.x}\t{c.y}\n" for c in cases)
        started = time.perf_counter()
        completed = subprocess.run(command, input=request, text=True,
                                   capture_output=True, timeout=180, check=False)
        summary = {"name": name, "mode": mode, "width": width, "epsilon": epsilon,
                   "cases": len(cases), "seconds": round(time.perf_counter() - started, 6)}
        self.batches.append(summary)
        assert completed.returncode == 0, (
            f"CLI exited {completed.returncode}: {completed.stderr[-2000:]}")
        lines = completed.stdout.splitlines()
        assert len(lines) == len(cases), f"expected {len(cases)} JSON lines, received {len(lines)}"
        for case, line in zip(cases, lines):
            self.checks += 1
            try:
                self.verify(json.loads(line), case, mode, epsilon, allow_budget_error)
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                if len(self.failures) < 50:
                    self.failures.append({"batch": name, "mode": mode, "case": case.name,
                                          "x": case.x, "y": case.y, "truth": case.truth,
                                          "reason": str(error), "output": line[:12000]})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "certalign")
    parser.add_argument("--output", type=Path, help="optionally write the JSON report here")
    args = parser.parse_args()
    started = time.perf_counter()
    audit = Audit(args.binary.resolve())
    exhaustive, stress = make_cases()
    try:
        for mode in MODES:
            audit.batch("exhaustive_binary_len_0_to_4", exhaustive, mode)
            audit.batch("random_and_adversarial", stress, mode)
        for mode in ("linear", "cost", "anchor-linear", "anchor-cost"):
            for width in (1, 3, 11):
                audit.batch("width_sweep", stress[:14] + stress[14::8], mode, width=width)
            for epsilon in ("0.1", "0.25", "1e-3", "1"):
                audit.batch("epsilon_guarantee", stress, mode, epsilon=epsilon)
            # A width-zero diagonal can cost 4 while the true distance is 2.
            # This request is strictly below the factor-two boundary despite
            # rounding to 1.0 in IEEE double precision.
            audit.batch("epsilon_decimal_boundary", [stress[14]], mode,
                        epsilon="0.999999999999999999")
        # This pair needs positive DP work in every mode. An exhausted budget
        # must not masquerade as an exact or approximate result.
        budget_case = [Case("budget_exhaustion", "ACGT" * 25, "TGCA" * 26,
                            oracle_distance("ACGT" * 25, "TGCA" * 26))]
        for mode in MODES:
            audit.batch("budget_exhaustion", budget_case, mode, max_cells=1,
                        allow_budget_error=True)
        assert audit.budget_errors == len(MODES), "some max-cells=1 runs ignored the budget"
    except (AssertionError, FileNotFoundError, subprocess.TimeoutExpired) as error:
        audit.failures.append({"reason": str(error), "kind": "batch_failure"})
    report = {"passed": not audit.failures,
              "oracle": "independent unrestricted unit-cost Levenshtein DP",
              "exhaustive_unique_pairs": len(exhaustive),
              "stress_unique_pairs": len(stress),
              "checked_outputs": audit.checks,
              "finite_trace_rounds": audit.finite_trace_rounds,
              "null_trace_rounds": audit.null_trace_rounds,
              "cigar_checks": audit.cigar_checks,
              "exact_results": audit.exact_results,
              "nonexact_approximate_results": audit.approximate_results,
              "budget_errors": audit.budget_errors,
              "seconds": round(time.perf_counter() - started, 6),
              "batches": audit.batches,
              "failures": audit.failures}
    rendered = json.dumps(report, indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(json.dumps({key: value for key, value in report.items() if key != "batches"}, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
