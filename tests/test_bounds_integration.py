#!/usr/bin/env python3
"""End-to-end bound audit with a separate full-DP oracle and v1 regression.

Exercises exact and rational-epsilon certification, CIGAR replay, finite
intermediate intervals, exhausted DP budgets, and monotonic strengthening.
The untouched v1 binary checks deterministic outputs, excluding timing fields.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from test_correctness import Audit, Case, MODES, make_cases

ROOT = Path(__file__).resolve().parents[1]
BOUNDS = ("length", "composition", "qgram")
V1_FIELDS = ("id", "mode", "n", "m", "distance", "lower_bound", "certified_exact",
             "certified_approx", "cigar", "rounds", "dp_cells", "certificate_edges",
             "probe_cells", "peak_cells", "trace", "error")


def run(binary: Path, cases: list[Case], mode: str, bound: str | None,
        width: int = 0, epsilon: str = "0", budget: int = 100_000_000) -> list[dict]:
    command = [str(binary), "--mode", mode, "--width", str(width), "--block", "64",
               "--seed-k", "11", "--epsilon", epsilon, "--max-cells", str(budget)]
    if bound is not None:
        command += ["--bound", bound]
    request = "".join(f"{case.name}\t{case.x}\t{case.y}\n" for case in cases)
    completed = subprocess.run(command, input=request, text=True, capture_output=True,
                               timeout=180, check=False)
    assert completed.returncode == 0, f"CLI failed: {completed.stderr[-2000:]}"
    results = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(results) == len(cases), f"expected {len(cases)} outputs, got {len(results)}"
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "certalign")
    parser.add_argument("--baseline", type=Path, default=ROOT / "reference_v1" / "build" / "certalign",
                        help="untouched v1 executable (default: reference_v1/build/certalign)")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    binary = args.binary.resolve()
    audit = Audit(binary)
    started = time.perf_counter()
    exhaustive, stress = make_cases()
    defaults_compared = baseline_compared = strengthening_comparisons = improved_rounds = 0
    budget_expectations = 0
    errors: list[dict] = []

    def verify_batch(name: str, cases: list[Case], mode: str, bound: str,
                     width: int = 0, epsilon: str = "0", budget: int = 100_000_000,
                     allow_error: bool = False) -> list[dict]:
        results = run(binary, cases, mode, bound, width, epsilon, budget)
        for case, result in zip(cases, results):
            audit.checks += 1
            try:
                audit.verify(result, case, mode, epsilon, allow_error)
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                if len(errors) < 20:
                    errors.append({"batch": name, "case": case.name, "mode": mode,
                                   "bound": bound, "reason": str(error), "output": result})
        return results

    try:
        for mode in MODES:
            cases = exhaustive + stress
            outputs = {bound: verify_batch("exact", cases, mode, bound) for bound in BOUNDS}
            for weak, strong in (("length", "composition"), ("composition", "qgram")):
                for case, a, b in zip(cases, outputs[weak], outputs[strong]):
                    assert b["dp_cells"] <= a["dp_cells"], f"stronger bound did more DP: {case.name}"
                    assert len(b["trace"]) <= len(a["trace"]), f"stronger bound did more rounds: {case.name}"
                    for old, new in zip(a["trace"], b["trace"]):
                        assert (new["width"], new["cells"], new["upper"]) == (
                            old["width"], old["cells"], old["upper"]), "predictor or DP changed"
                        if old["lower"] is not None and new["lower"] is not None:
                            assert new["lower"] >= old["lower"], "strengthened lower bound decreased"
                            improved_rounds += new["lower"] > old["lower"]
                        strengthening_comparisons += 1
            default_results = run(binary, cases, mode, None)
            for default, explicit in zip(default_results, outputs["length"]):
                assert {k: default.get(k) for k in V1_FIELDS} == {
                    k: explicit.get(k) for k in V1_FIELDS}, "default differs from explicit length bound"
                defaults_compared += 1
            if args.baseline:
                old_results = run(args.baseline.resolve(), cases, mode, None)
                for old, new in zip(old_results, outputs["length"]):
                    assert {k: old.get(k) for k in V1_FIELDS} == {
                        k: new.get(k) for k in V1_FIELDS}, "length mode differs from untouched v1"
                    baseline_compared += 1
        selected = stress[:22] + stress[22::9]
        for bound in BOUNDS:
            for mode in ("linear", "cost", "anchor-cost"):
                for width in (1, 7):
                    verify_batch("width", selected, mode, bound, width=width)
                for epsilon in ("0.1", "0.25", "1", "0.999999999999999999"):
                    verify_batch("epsilon", selected, mode, bound, epsilon=epsilon)
            # There is a finite first round, but it cannot certify the wrong
            # width-zero score 4 because the unrestricted optimum is 2.
            case = Case("budget_after_finite_round", "ACAC", "CACA", 2)
            for budget in (0, 1, 5):
                result = verify_batch("budget", [case], "linear", bound, budget=budget,
                                      allow_error=True)[0]
                assert result.get("error"), "exhausted budget was not reported"
                assert len(result["trace"]) == (1 if budget == 5 else 0)
                budget_expectations += 1
            for budget in (24, 25):
                result = verify_batch("full_budget_boundary", [case], "full", bound,
                                      budget=budget, allow_error=budget == 24)[0]
                assert bool(result.get("error")) == (budget == 24)
                budget_expectations += 1
            empty = Case("empty_budget_boundary", "", "", 0)
            for budget in (0, 1):
                result = verify_batch("empty_budget_boundary", [empty], "linear", bound,
                                      budget=budget, allow_error=budget == 0)[0]
                assert bool(result.get("error")) == (budget == 0)
                budget_expectations += 1
        invalid = subprocess.run([str(binary), "--bound", "invalid"], input="", text=True,
                                 capture_output=True, timeout=10, check=False)
        assert invalid.returncode != 0 and not invalid.stdout, "invalid bound was accepted"
        assert improved_rounds > 0, "bounds did not strengthen any tested certificate"
    except (AssertionError, FileNotFoundError, subprocess.TimeoutExpired) as error:
        errors.append({"reason": str(error), "kind": "batch_failure"})
    report = {"passed": not errors, "oracle": "independent unrestricted unit-cost Levenshtein DP",
              "checked_outputs": audit.checks, "finite_trace_rounds": audit.finite_trace_rounds,
              "null_trace_rounds": audit.null_trace_rounds, "cigar_checks": audit.cigar_checks,
              "exact_results": audit.exact_results,
              "nonexact_approximate_results": audit.approximate_results,
              "budget_errors": audit.budget_errors, "budget_expectations": budget_expectations,
              "default_length_comparisons": defaults_compared,
              "v1_deterministic_comparisons": baseline_compared,
              "v1_regression_tested": args.baseline is not None,
              "strengthening_round_comparisons": strengthening_comparisons,
              "strictly_improved_round_bounds": improved_rounds,
              "seconds": round(time.perf_counter() - started, 6), "failures": errors}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
