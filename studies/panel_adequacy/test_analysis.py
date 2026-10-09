#!/usr/bin/env python3
"""Independent toy checks for leakage, calibration, ties, and clustered sampling."""
import math
import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from analyze import (Predictors, Truth, FrozenLogistic, feature_matrix, split_ids,
                     stable_hash, matched_accept, calibration_cutoff,
                     threshold_accept, metric, bootstrap_weights, scores, read_cosigt)


def case(sample, locus, q=10, n=100, u=2, w=1):
    return Predictors(sample, locus, q, n, u, w, "*")


class AnalysisTests(unittest.TestCase):
    def test_identity_split_reproducible_disjoint(self):
        ids = [str(i) for i in range(40)]
        result = split_ids(ids, "donor")
        self.assertEqual(result, split_ids(list(reversed(ids)) + ids, "donor"))
        self.assertEqual([list(result.values()).count(k) for k in ("development", "calibration", "test")], [20, 10, 10])

    def test_ties_cutoff_is_same_as_calibration_rank(self):
        rows = [case(str(i), "L") for i in range(20)]
        risk = np.zeros(20)
        ranked = matched_accept(risk, rows, .95)
        self.assertEqual(int(ranked.sum()), 19)
        self.assertTrue(np.array_equal(ranked, threshold_accept(risk, rows, calibration_cutoff(risk, rows))))
        excluded = np.where(~ranked)[0][0]
        self.assertEqual(rows[excluded].tie, max(r.tie for r in rows))

    def test_cutoff_transfer_is_not_test_quantile(self):
        cal = [case(str(i), "C") for i in range(20)]
        test = [case(str(i), "T") for i in range(10)]
        cutoff = calibration_cutoff(np.arange(20.), cal)
        self.assertEqual(int(threshold_accept(np.arange(100., 110.), test, cutoff).sum()), 0)
        self.assertEqual(int(matched_accept(np.arange(100., 110.), test, .95).sum()), 10)

    def test_feature_whitelist_and_missing(self):
        self.assertEqual(set(Predictors.__dataclass_fields__),
                         {"sample", "locus", "quality", "total_reads", "unexpl_reads", "weight_dist", "warnings"})
        x = feature_matrix([case("A", "L", n=0), case("B", "L", q=-1)])
        self.assertTrue(np.isnan(x[0, 4]))
        self.assertTrue(np.isnan(x[1, 0]))
        self.assertAlmostEqual(x[1, 2], math.log1p(2))

    def test_model_preprocessing_never_refits_on_prediction(self):
        train = [case(str(i), "L", q=i, w=10-i) for i in range(10)]
        labels = np.asarray([i < 5 for i in range(10)])
        model = FrozenLogistic().fit(train, labels)
        means, medians, coefficients = model.scaler.mean_.copy(), model.medians.copy(), model.model.coef_.copy()
        a = model.predict([case("X", "L", q=1)])[0]
        b = model.predict([case("X", "L", q=1), case("Y", "L", q=1e20)])[0]
        self.assertAlmostEqual(a, b, places=14)
        self.assertTrue(np.array_equal(model.scaler.mean_, means))
        self.assertTrue(np.array_equal(model.medians, medians))
        self.assertTrue(np.array_equal(model.model.coef_, coefficients))

    def test_error_metric_counts_errors_only_when_retained(self):
        labels = np.asarray([True, True, False])
        keep = np.asarray([True, False, True])
        m = metric(labels, keep, np.asarray([2, 5, 3]))
        self.assertEqual(m["accepted_wrong"], 2)
        self.assertEqual(m["accepted"], 5)
        self.assertEqual(m["error_per_eligible"], .2)
        self.assertEqual(m["error_per_accepted"], .4)

    def test_cluster_bootstrap_keeps_overlap_loci_together(self):
        rows = [case(s, l) for s in ("A", "B", "C") for l in ("L1", "L2", "L3")]
        groups = {"L1": "G1", "L2": "G1", "L3": "G2"}
        ws = list(bootstrap_weights(rows, groups, 25, 99))
        self.assertTrue(any(not np.array_equal(w, np.ones(9)) for w in ws))
        for w in ws:
            for i in (0, 3, 6):
                self.assertEqual(w[i], w[i+1])
        self.assertTrue(all(np.array_equal(a,b) for a,b in zip(ws, bootstrap_weights(rows, groups, 25, 99))))

    def test_native_risk_directions(self):
        train = [case(str(i), "L", q=i, w=10-i) for i in range(10)]
        model = FrozenLogistic().fit(train, np.asarray([i < 5 for i in range(10)]))
        rows = [case("A", "L", q=20, n=100, u=1, w=1), case("B", "L", q=2, n=10, u=2, w=3)]
        risks = scores(rows, model)
        for name in ("quality", "total_reads", "unexpl_reads", "weight_dist", "unexpl_fraction"):
            self.assertLess(risks[name][0], risks[name][1])

    def test_truth_validity_preserves_infinite_qv(self):
        self.assertTrue(Truth("A", "L", 0, 100, 0, np.inf, 0, 100, 0, np.inf).scorable)
        for size, div in [(0, 0), (100, -1), (100, np.nan), (np.inf, 0)]:
            self.assertFalse(Truth("A", "L", 0, size, div, np.inf, 0, 100, 0, np.inf).scorable)

    def test_one_class_is_unavailable_without_changing_other_scores(self):
        rows = [case(str(i), "L", q=i) for i in range(5)]
        model = FrozenLogistic().fit(rows, np.zeros(5, dtype=bool))
        self.assertFalse(model.metadata()["available"])
        self.assertNotIn("logistic", scores(rows, model))
        self.assertIn("quality", scores(rows, model))

    def test_panel_bootstrap_is_paired_with_missing_cluster(self):
        rows = [case(s, l) for s in ("A", "B", "C") for l in ("L1", "L2")]
        groups = {"L1": "G1", "L2": "G2"}
        all_weights = list(bootstrap_weights(rows, groups, 20, 99, ["A", "B", "C"], ["G1", "G2"]))
        subset_weights = list(bootstrap_weights(rows[:2], groups, 20, 99, ["A", "B", "C"], ["G1", "G2"]))
        self.assertTrue(all(np.array_equal(a[:2], b) for a, b in zip(all_weights, subset_weights)))

    def test_cosigt_deduplication_refuses_conflicting_rows(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "toy.tsv.gz"
            with gzip.open(path, "wt") as fh:
                fh.write("sample\tregion\tscore\nA\tL\t1\nA\tL\t1\n")
            rows, meta = read_cosigt(path)
            self.assertEqual(len(rows), 1)
            self.assertEqual(meta["byte_identical_duplicates_removed"], 1)
            with gzip.open(path, "wt") as fh:
                fh.write("sample\tregion\tscore\nA\tL\t1\nA\tL\t2\n")
            with self.assertRaises(ValueError):
                read_cosigt(path)


if __name__ == "__main__":
    unittest.main()
