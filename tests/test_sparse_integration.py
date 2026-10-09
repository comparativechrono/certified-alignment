#!/usr/bin/env python3
"""Independent sparse-certificate integration and unchanged-v3 regression audit.

Reuses an independent unrestricted edit-distance oracle and CIGAR replayer.
The capped variant must remain at least as effective as its retained qgram
bound, while uncapped exact certification must stop at the dense bound's round.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

from test_correctness import Audit, Case, MODES, make_cases
from test_ordered_integration import DETERMINISTIC, V2_BOUND_FIELDS, run

ROOT = Path(__file__).resolve().parents[1]
V3_BOUND_FIELDS = V2_BOUND_FIELDS + ("ordered_cells", "ordered_bytes", "ordered_calls",
                                    "ordered_fallback")
EXACT_FIELDS = ("id", "mode", "n", "m", "distance", "certified_exact",
                "certified_approx", "cigar", "rounds", "dp_cells", "peak_cells", "error")
UNCAPPED = ("--sparse-max-work", "18446744073709551615")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "build" / "certalign")
    parser.add_argument("--baseline", type=Path, default=ROOT / "reference_v3" / "build" / "certalign")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    binary = args.binary.resolve()
    audit = Audit(binary)
    exhaustive, stress = make_cases()
    all_cases = exhaustive + stress
    selected = exhaustive[::31] + stress[:22] + stress[22::9]
    baseline_compared = dense_compared = strengthening_rounds = 0
    strict_improvements = fallback_compared = budget_expectations = active_outputs = 0
    cap_outputs = 0
    errors: list[dict] = []

    def verify(name: str, cases: list[Case], mode: str, bound: str, **kwargs) -> list[dict]:
        nonlocal active_outputs, cap_outputs
        allow_error = kwargs.pop("allow_error", False)
        results = run(binary, cases, mode, bound, **kwargs)
        extra = kwargs.get("extra", ())
        cap = int(extra[extra.index("--sparse-max-work") + 1]) if "--sparse-max-work" in extra else 4_000_000
        memory = int(extra[extra.index("--bound-memory") + 1]) if "--bound-memory" in extra else 64 * 1024 * 1024
        for case, result in zip(cases, results):
            audit.checks += 1
            active_outputs += "+sparse" in result.get("bound_active", "")
            try:
                audit.verify(result, case, mode, kwargs.get("epsilon", "0"), allow_error)
                assert result.get("sparse_expansions", 0) <= cap, "sparse work cap exceeded"
                assert result.get("bound_bytes", 0) <= memory, "shared bound memory exceeded"
                cap_outputs += 1
            except (AssertionError, KeyError, TypeError, ValueError) as error:
                if len(errors) < 20:
                    errors.append({"batch": name, "case": case.name, "mode": mode,
                                   "bound": bound, "reason": str(error), "output": result})
        return results

    def compare_stronger(cases: list[Case], weak: list[dict], strong: list[dict]) -> None:
        nonlocal strengthening_rounds, strict_improvements
        for case, a, b in zip(cases, weak, strong):
            assert b["dp_cells"] <= a["dp_cells"], f"sparse used more DP than qgram: {case.name}"
            assert len(b["trace"]) <= len(a["trace"]), f"sparse used more rounds: {case.name}"
            for old, new in zip(a["trace"], b["trace"]):
                assert (old["width"], old["cells"], old["upper"]) == (
                    new["width"], new["cells"], new["upper"]), "predictor or restricted DP changed"
                if old["lower"] is not None and new["lower"] is not None:
                    assert new["lower"] >= old["lower"], f"retained qgram bound decreased: {case.name}"
                    strict_improvements += new["lower"] > old["lower"]
                strengthening_rounds += 1

    def compare_dense(cases: list[Case], dense: list[dict], sparse: list[dict]) -> None:
        nonlocal dense_compared
        for case, a, b in zip(cases, dense, sparse):
            assert {k: a.get(k) for k in EXACT_FIELDS} == {k: b.get(k) for k in EXACT_FIELDS}, (
                f"uncapped exact sparse/dense mismatch: {case.name}")
            for old, new in zip(a["trace"], b["trace"]):
                assert (old["width"], old["cells"], old["upper"]) == (
                    new["width"], new["cells"], new["upper"]), "sparse changed DP geometry"
                if old["lower"] is not None and new["lower"] is not None:
                    assert new["lower"] <= old["lower"], "sparse exceeded dense relaxed bound"
            dense_compared += 1

    try:
        for mode in MODES:
            q_results = dense_results = None
            for bound in ("length", "composition", "qgram", "ordered"):
                new = verify("original_bounds", all_cases, mode, bound)
                old = run(args.baseline.resolve(), all_cases, mode, bound)
                for a, b in zip(old, new):
                    assert {k: a.get(k) for k in DETERMINISTIC + V3_BOUND_FIELDS} == {
                        k: b.get(k) for k in DETERMINISTIC + V3_BOUND_FIELDS}, (
                        f"v3 deterministic output changed: {mode}/{bound}/{a.get('id')}")
                    baseline_compared += 1
                if bound == "qgram": q_results = new
                if bound == "ordered": dense_results = new
            default = verify("sparse_default", all_cases, mode, "sparse")
            compare_stronger(all_cases, q_results, default)
            uncapped = verify("sparse_uncapped", all_cases, mode, "sparse", extra=UNCAPPED)
            compare_dense(all_cases, dense_results, uncapped)
            for k in range(1, 8):
                old = run(binary, selected, mode, "ordered", extra=("--order-block", str(k)))
                new = verify("sparse_block_sweep", selected, mode, "sparse",
                             extra=UNCAPPED + ("--order-block", str(k)))
                compare_dense(selected, old, new)

        for mode in ("linear", "cost", "anchor-cost"):
            for k in (1, 3, 8):
                for epsilon in ("0.1", "0.25", "1", "0.999999999999999999"):
                    old = run(binary, selected, mode, "qgram", epsilon=epsilon)
                    new = verify("epsilon", selected, mode, "sparse", epsilon=epsilon,
                                 extra=("--order-block", str(k)))
                    compare_stronger(selected, old, new)
            for width in (1, 7):
                old = run(binary, selected, mode, "ordered", width=width)
                new = verify("width", selected, mode, "sparse", width=width, extra=UNCAPPED)
                compare_dense(selected, old, new)
            for work in (0, 1, 7, 64):
                old = run(binary, selected, mode, "qgram")
                new = verify("work_cap", selected, mode, "sparse",
                             extra=("--sparse-max-work", str(work)))
                compare_stronger(selected, old, new)
                if not work:
                    for a, b in zip(old, new):
                        assert {k: a.get(k) for k in DETERMINISTIC} == {
                            k: b.get(k) for k in DETERMINISTIC}, "zero work changed retained qgram behavior"
                        assert b.get("sparse_expansions", 0) == 0
                        fallback_compared += 1
            for memory in (0, 1, 4096, 65536):
                old = run(binary, selected, mode, "qgram", extra=("--bound-memory", str(memory)))
                new = verify("memory_cap", selected, mode, "sparse",
                             extra=("--bound-memory", str(memory)))
                compare_stronger(selected, old, new)
                if not memory:
                    for a, b in zip(old, new):
                        assert {k: a.get(k) for k in DETERMINISTIC} == {
                            k: b.get(k) for k in DETERMINISTIC}, "zero memory changed length-only behavior"
                        fallback_compared += 1
        case = Case("budget_after_finite_round", "ACAC", "CACA", 2)
        for budget in (0, 1, 5):
            result = verify("budget", [case], "linear", "sparse", budget=budget,
                            allow_error=True)[0]
            assert result.get("error"), "exhausted budget not reported"
            assert len(result["trace"]) == (1 if budget == 5 else 0)
            budget_expectations += 1
        for budget in (24, 25):
            result = verify("full_budget", [case], "full", "sparse", budget=budget,
                            allow_error=budget == 24)[0]
            assert bool(result.get("error")) == (budget == 24)
            budget_expectations += 1
        for extra in (("--order-block", "0"), ("--order-block", "9"),
                      ("--sparse-max-work", "not-a-number"), ("--sparse-max-work", "-1")):
            invalid = subprocess.run([str(binary), "--bound", "sparse", *extra], input="",
                                     text=True, capture_output=True, timeout=10, check=False)
            assert invalid.returncode != 0 and not invalid.stdout, f"invalid flag accepted: {extra}"
        assert strict_improvements > 0, "sparse strengthened no certificate"
        assert active_outputs > 0, "sparse cache never activated"
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
        "v3_deterministic_comparisons": baseline_compared,
        "uncapped_dense_exact_comparisons": dense_compared,
        "strengthening_round_comparisons": strengthening_rounds,
        "strictly_improved_round_bounds": strict_improvements,
        "fallback_comparisons": fallback_compared,
        "outputs_with_sparse_index": active_outputs,
        "outputs_checked_for_caps": cap_outputs,
        "allocation_failure_injection_tested": False,
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
