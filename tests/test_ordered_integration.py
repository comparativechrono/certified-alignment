#!/usr/bin/env python3
"""Independent end-to-end ordered-bound and unchanged-v2 regression audit.

Reuses the separately implemented unrestricted edit-distance oracle and CIGAR
replayer. The v2 executable checks every original bound/predictor combination;
the ordered bound is compared with qgram on identical predicted bands.
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
DETERMINISTIC = ("id", "mode", "n", "m", "distance", "lower_bound", "certified_exact",
                 "certified_approx", "cigar", "rounds", "dp_cells", "certificate_edges",
                 "probe_cells", "peak_cells", "trace", "error")
V2_BOUND_FIELDS = ("bound", "bound_bytes", "bound_calls", "composition_calls", "q2_calls",
                   "q4_calls", "bound_cursor_updates", "bound_active", "bound_fallback")


def run(binary: Path, cases: list[Case], mode: str, bound: str,
        *, width: int = 0, epsilon: str = "0", budget: int = 100_000_000,
        extra: tuple[str, ...] = ()) -> list[dict]:
    command = [str(binary), "--mode", mode, "--width", str(width), "--block", "64",
               "--seed-k", "11", "--epsilon", epsilon, "--max-cells", str(budget),
               "--bound", bound, *extra]
    completed = subprocess.run(command,
        input="".join(f"{case.name}\t{case.x}\t{case.y}\n" for case in cases),
        text=True, capture_output=True, timeout=180, check=False)
    assert completed.returncode == 0, f"CLI failed: {completed.stderr[-2000:]}"
    results = [json.loads(line) for line in completed.stdout.splitlines()]
    assert len(results) == len(cases), f"expected {len(cases)} outputs, got {len(results)}"
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "certalign")
    parser.add_argument("--baseline", type=Path, default=ROOT / "reference_v2" / "build" / "certalign")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    binary = args.binary.resolve()
    audit = Audit(binary)
    exhaustive, stress = make_cases()
    all_cases = exhaustive + stress
    selected = exhaustive[::31] + stress[:22] + stress[22::9]
    baseline_compared = strengthening_rounds = strict_improvements = 0
    fallback_compared = budget_expectations = 0
    active_outputs = 0
    errors: list[dict] = []

    def verify(name: str, cases: list[Case], mode: str, bound: str, **kwargs) -> list[dict]:
        nonlocal active_outputs
        allow_error = kwargs.pop("allow_error", False)
        results = run(binary, cases, mode, bound, **kwargs)
        for case, result in zip(cases, results):
            audit.checks += 1
            active_outputs += "+ordered" in result.get("bound_active", "")
            try:
                audit.verify(result, case, mode, kwargs.get("epsilon", "0"), allow_error)
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                if len(errors) < 20:
                    errors.append({"batch": name, "case": case.name, "mode": mode,
                                   "bound": bound, "reason": str(error), "output": result})
        return results

    def compare_stronger(cases: list[Case], old: list[dict], new: list[dict]) -> None:
        nonlocal strengthening_rounds, strict_improvements
        for case, a, b in zip(cases, old, new):
            assert b["dp_cells"] <= a["dp_cells"], f"ordered used more DP cells: {case.name}"
            assert len(b["trace"]) <= len(a["trace"]), f"ordered used more rounds: {case.name}"
            for weak, strong in zip(a["trace"], b["trace"]):
                assert (weak["width"], weak["cells"], weak["upper"]) == (
                    strong["width"], strong["cells"], strong["upper"]), "predictor or restricted DP changed"
                if weak["lower"] is not None and strong["lower"] is not None:
                    assert strong["lower"] >= weak["lower"], f"ordered lower bound decreased: {case.name}"
                    strict_improvements += strong["lower"] > weak["lower"]
                strengthening_rounds += 1

    try:
        for mode in MODES:
            q_results = None
            for bound in ("length", "composition", "qgram"):
                new = verify("original_bounds", all_cases, mode, bound)
                old = run(args.baseline.resolve(), all_cases, mode, bound)
                for a, b in zip(old, new):
                    assert {k: a.get(k) for k in DETERMINISTIC + V2_BOUND_FIELDS} == {
                        k: b.get(k) for k in DETERMINISTIC + V2_BOUND_FIELDS}, (
                        f"v2 deterministic output changed: {mode}/{bound}/{a.get('id')}")
                    baseline_compared += 1
                if bound == "qgram":
                    q_results = new
            ordered = verify("ordered_default", all_cases, mode, "ordered")
            compare_stronger(all_cases, q_results, ordered)
            for k in range(1, 8):
                old = run(binary, selected, mode, "qgram")
                new = verify("ordered_block_sweep", selected, mode, "ordered",
                             extra=("--order-block", str(k)))
                compare_stronger(selected, old, new)
        for mode in ("linear", "cost", "anchor-cost"):
            for k in (1, 3, 8):
                for epsilon in ("0.1", "0.25", "1", "0.999999999999999999"):
                    old = run(binary, selected, mode, "qgram", epsilon=epsilon)
                    new = verify("epsilon", selected, mode, "ordered", epsilon=epsilon,
                                 extra=("--order-block", str(k)))
                    compare_stronger(selected, old, new)
            for width in (1, 7):
                old = run(binary, selected, mode, "qgram", width=width)
                new = verify("width", selected, mode, "ordered", width=width)
                compare_stronger(selected, old, new)
            for extra, reference in ((("--order-max-cells", "0"), "qgram"),
                                     (("--order-max-cells", "1"), "qgram"),
                                     (("--bound-memory", "0"), "length")):
                old = run(binary, selected, mode, reference)
                new = verify("fallback", selected, mode, "ordered", extra=extra)
                for a, b in zip(old, new):
                    assert {k: a.get(k) for k in DETERMINISTIC} == {
                        k: b.get(k) for k in DETERMINISTIC}, "fallback changed retained bound behavior"
                    assert not b.get("ordered_cells"), "fallback built ordered table"
                    fallback_compared += 1
        case = Case("budget_after_finite_round", "ACAC", "CACA", 2)
        for budget in (0, 1, 5):
            result = verify("budget", [case], "linear", "ordered", budget=budget,
                            allow_error=True)[0]
            assert result.get("error"), "exhausted budget was not reported"
            assert len(result["trace"]) == (1 if budget == 5 else 0)
            budget_expectations += 1
        for budget in (24, 25):
            result = verify("full_budget", [case], "full", "ordered", budget=budget,
                            allow_error=budget == 24)[0]
            assert bool(result.get("error")) == (budget == 24)
            budget_expectations += 1
        for extra in (("--order-block", "0"), ("--order-block", "9"),
                      ("--order-max-cells", "not-a-number")):
            invalid = subprocess.run([str(binary), "--bound", "ordered", *extra], input="",
                                     text=True, capture_output=True, timeout=10, check=False)
            assert invalid.returncode != 0 and not invalid.stdout, f"invalid flag accepted: {extra}"
        assert strict_improvements > 0, "ordered bound strengthened no tested certificates"
        assert active_outputs > 0, "ordered index never activated"
    except (AssertionError, FileNotFoundError, subprocess.TimeoutExpired) as error:
        errors.append({"kind": "batch_failure", "reason": str(error)})
    report = {
        "passed": not errors,
        "oracle": "independent unrestricted unit-cost Levenshtein DP and independent CIGAR replay",
        "checked_outputs": audit.checks,
        "finite_trace_rounds": audit.finite_trace_rounds,
        "null_trace_rounds": audit.null_trace_rounds,
        "cigar_checks": audit.cigar_checks,
        "exact_results": audit.exact_results,
        "nonexact_approximate_results": audit.approximate_results,
        "budget_errors": audit.budget_errors,
        "budget_expectations": budget_expectations,
        "v2_deterministic_comparisons": baseline_compared,
        "strengthening_round_comparisons": strengthening_rounds,
        "strictly_improved_round_bounds": strict_improvements,
        "fallback_comparisons": fallback_compared,
        "outputs_with_ordered_index": active_outputs,
        "seconds": round(time.perf_counter() - started, 6),
        "failures": errors,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
