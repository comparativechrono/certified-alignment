#!/usr/bin/env python3
"""Frozen retrospective native-confidence baselines; no aligner is rerun.

Inputs, fitting, and decisions are separated from truth evaluation. The fit uses
only development LOO labels. Native-baseline selection uses calibration labels.
Neither test labels nor full-panel labels enter fitting, ranking, or cutoffs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import scipy
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


SEED = 20261009
BOOTSTRAPS = 4999
RETENTIONS = (0.8, 0.9, 0.95, 0.99, 1.0)
ERROR_THRESHOLDS = (0.005, 0.001, 0.01)
NATIVE_METHODS = ("quality", "weight_dist", "unexpl_fraction", "unexpl_reads", "total_reads")
FEATURE_NAMES = ("log1p_quality", "log1p_weight_dist", "log1p_unexpl_reads",
                 "log1p_total_reads", "unexpl_fraction")
TRUTH_NAMES = ("edit", "size", "div", "qv", "avail_edit", "avail_size", "avail_div", "avail_qv")


def stable_hash(kind: str, value: str) -> str:
    return hashlib.sha256(f"panel-native-v1|{kind}|{value}".encode()).hexdigest()


def split_ids(ids: list[str], kind: str) -> dict[str, str]:
    ordered = sorted(set(ids), key=lambda value: stable_hash(kind, value))
    n_dev, n_cal = math.floor(len(ordered) * .5), math.floor(len(ordered) * .25)
    return {value: ("development" if i < n_dev else
                    "calibration" if i < n_dev + n_cal else "test")
            for i, value in enumerate(ordered)}


def number(value: str | None) -> float:
    try:
        return float(value) if value is not None else float("nan")
    except ValueError:
        return float("nan")


@dataclass(frozen=True)
class Predictors:
    sample: str
    locus: str
    quality: float
    total_reads: float
    unexpl_reads: float
    weight_dist: float
    warnings: str

    @property
    def key(self) -> tuple[str, str]:
        return self.sample, self.locus

    @property
    def tie(self) -> str:
        return stable_hash("tie", self.sample + "|" + self.locus)

    @property
    def unexpl_fraction(self) -> float:
        if not (np.isfinite(self.total_reads) and self.total_reads > 0
                and np.isfinite(self.unexpl_reads)):
            return float("nan")
        return self.unexpl_reads / self.total_reads


@dataclass(frozen=True)
class Truth:
    sample: str
    locus: str
    edit: float
    size: float
    div: float
    qv: float
    avail_edit: float
    avail_size: float
    avail_div: float
    avail_qv: float

    @property
    def key(self) -> tuple[str, str]:
        return self.sample, self.locus

    @property
    def scorable(self) -> bool:
        return bool(np.isfinite(self.div) and self.div >= 0 and np.isfinite(self.size) and self.size > 0)


def read_locityper(path: Path) -> tuple[dict, dict, dict]:
    predictors, truth = {}, {}
    comments, total_rows = [], 0
    with gzip.open(path, "rt") as fh:
        def noncomments():
            for line in fh:
                if line.startswith("#"):
                    comments.append(line.rstrip())
                else:
                    yield line
        reader = csv.DictReader(noncomments(), delimiter="\t")
        for row in reader:
            total_rows += 1
            if row["query_type"] != "gt":
                continue
            p = Predictors(row["sample"], row["locus"], number(row["quality"]),
                           number(row["total_reads"]), number(row["unexpl_reads"]),
                           number(row["weight_dist"]), row["warnings"])
            if any(np.isfinite(v) and v < 0 for v in (p.quality, p.total_reads, p.unexpl_reads, p.weight_dist)):
                raise ValueError(f"Negative native predictor in {path}: {p.key}")
            if p.key in predictors:
                raise ValueError(f"Duplicate genotype key in {path}: {p.key}")
            predictors[p.key] = p
            truth[p.key] = Truth(p.sample, p.locus, *(number(row.get(k)) for k in TRUTH_NAMES))
    return predictors, truth, {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                              "rows": total_rows, "gt_rows": len(predictors), "comments": comments}


def feature_matrix(rows: list[Predictors]) -> np.ndarray:
    values = []
    for p in rows:
        raw = (p.quality, p.weight_dist, p.unexpl_reads, p.total_reads)
        transformed = [math.log1p(v) if np.isfinite(v) and v >= 0 else np.nan for v in raw]
        ratio = p.unexpl_fraction
        transformed.append(ratio if np.isfinite(ratio) and ratio >= 0 else np.nan)
        values.append(transformed)
    return np.asarray(values, dtype=float).reshape((-1, len(FEATURE_NAMES)))


class FrozenLogistic:
    def fit(self, predictors: list[Predictors], labels: np.ndarray):
        x = feature_matrix(predictors)
        self.medians = np.asarray([np.median(col[np.isfinite(col)]) if np.isfinite(col).any() else 0.
                                   for col in x.T])
        z = self.impute(x)
        self.scaler = StandardScaler().fit(z)
        if len(np.unique(labels)) < 2:
            self.model = None
            self.unavailable_reason = "Development labels have fewer than two classes"
            return self
        self.model = LogisticRegression(C=1., penalty="l2", max_iter=2000,
                                        class_weight=None, random_state=SEED, solver="lbfgs")
        self.model.fit(self.scaler.transform(z), labels)
        if np.any(self.model.n_iter_ >= 2000):
            raise RuntimeError("Frozen logistic regression did not converge")
        return self

    def impute(self, x):
        missing = ~np.isfinite(x)
        return np.concatenate([np.where(missing, self.medians, x), missing.astype(float)], axis=1)

    def predict(self, predictors: list[Predictors]) -> np.ndarray:
        if self.model is None:
            raise ValueError(self.unavailable_reason)
        return self.model.predict_proba(self.scaler.transform(self.impute(feature_matrix(predictors))))[:, 1]

    def metadata(self):
        names = list(FEATURE_NAMES) + [name + "_missing" for name in FEATURE_NAMES]
        if self.model is None:
            return {"available": False, "unavailable_reason": self.unavailable_reason, "feature_names": names}
        return {"available": True, "feature_names": names, "imputation_medians": self.medians.tolist(),
                "scaler_mean": self.scaler.mean_.tolist(), "scaler_scale": self.scaler.scale_.tolist(),
                "coefficients_standardized": self.model.coef_[0].tolist(),
                "intercept": self.model.intercept_[0], "iterations": int(self.model.n_iter_[0]),
                "C": 1., "penalty": "l2", "solver": "lbfgs", "max_iter": 2000,
                "class_weight": None, "random_state": SEED}


def scores(predictors: list[Predictors], model: FrozenLogistic) -> dict[str, np.ndarray]:
    def finite_risk(value):
        return value if np.isfinite(value) else float("inf")
    result = {
        "quality": np.asarray([finite_risk(-p.quality) for p in predictors]),
        "weight_dist": np.asarray([finite_risk(p.weight_dist) for p in predictors]),
        "unexpl_fraction": np.asarray([finite_risk(p.unexpl_fraction) for p in predictors]),
        "unexpl_reads": np.asarray([finite_risk(p.unexpl_reads) for p in predictors]),
        "total_reads": np.asarray([finite_risk(-p.total_reads) for p in predictors]),
        "hash_control": np.zeros(len(predictors)),
        "no_filter": np.zeros(len(predictors)),
    }
    if model.model is not None:
        result["logistic"] = model.predict(predictors)
    return result


def rank_order(risk: np.ndarray, predictors: list[Predictors]) -> np.ndarray:
    return np.asarray(sorted(range(len(predictors)), key=lambda i: (risk[i], predictors[i].tie)), dtype=int)


def matched_accept(risk, predictors, retention):
    n = math.ceil(retention * len(predictors))
    accept = np.zeros(len(predictors), dtype=bool)
    accept[rank_order(risk, predictors)[:n]] = True
    return accept


def calibration_cutoff(risk, predictors, retention=.95):
    if not predictors:
        raise ValueError("Empty calibration data")
    i = rank_order(risk, predictors)[math.ceil(retention * len(predictors)) - 1]
    return float(risk[i]), predictors[i].tie


def threshold_accept(risk, predictors, cutoff):
    return np.asarray([(float(value), p.tie) <= cutoff for value, p in zip(risk, predictors)])


def metric(labels, accept, weights=None):
    w = np.ones(len(labels)) if weights is None else np.asarray(weights)
    n, retained, wrong = float(w.sum()), float(w[accept].sum()), float(w[accept & labels].sum())
    return {"eligible": n, "accepted": retained, "coverage": retained / n if n else None,
            "accepted_wrong": wrong, "error_per_eligible": wrong / n if n else None,
            "error_per_accepted": wrong / retained if retained else None,
            "total_wrong": float(w[labels].sum()),
            "fraction_errors_removed": 1 - wrong / float(w[labels].sum()) if w[labels].sum() else None}


def quantiles(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    return {"n": len(values), "finite": len(finite), "nonfinite": len(values) - len(finite),
            "mean": float(finite.mean()) if len(finite) else None,
            "quantiles": dict(zip(("0", ".5", ".9", ".95", ".99", "1"),
                                  np.quantile(finite, [0, .5, .9, .95, .99, 1]).tolist())) if len(finite) else {}}


def describe(truth):
    rows = list(truth.values())
    div = np.asarray([r.div for r in rows])
    ok = np.asarray([r.scorable for r in rows])
    available = np.asarray([r.avail_div for r in rows])
    valid_pair = ok & np.isfinite(available)
    return {"gt_rows": len(rows), "truth_missing_excluded": int((~ok).sum()),
            "divergence": quantiles(div[ok]),
            "wrong_counts": {str(t): int((div[ok] > t).sum()) for t in ERROR_THRESHOLDS},
            "closest_available_divergence": quantiles(available),
            "called_minus_available_divergence": quantiles(div[valid_pair] - available[valid_pair]),
            "called_below_available_count": int((div[valid_pair] < available[valid_pair]).sum()),
            "available_is_not_a_strict_floor": True}


def bootstrap_weights(predictors, groups, replicates=BOOTSTRAPS, seed=SEED,
                      donor_universe=None, group_universe=None):
    """Independent cluster resampling along donor and overlap-group dimensions."""
    donor_ids = sorted(donor_universe if donor_universe is not None else {p.sample for p in predictors})
    group_ids = sorted(group_universe if group_universe is not None else {groups[p.locus] for p in predictors})
    donor_index = {s: i for i, s in enumerate(donor_ids)}
    group_index = {s: i for i, s in enumerate(group_ids)}
    di = np.asarray([donor_index[p.sample] for p in predictors])
    gi = np.asarray([group_index[groups[p.locus]] for p in predictors])
    rng = np.random.default_rng(seed)
    for _ in range(replicates):
        donor_counts = np.bincount(rng.integers(len(donor_ids), size=len(donor_ids)), minlength=len(donor_ids))
        group_counts = np.bincount(rng.integers(len(group_ids), size=len(group_ids)), minlength=len(group_ids))
        yield donor_counts[di] * group_counts[gi]


def percentile_interval(values):
    finite = np.asarray([v for v in values if v is not None and np.isfinite(v)])
    return {"low": float(np.quantile(finite, .025)) if len(finite) else None,
            "high": float(np.quantile(finite, .975)) if len(finite) else None,
            "valid_replicates": len(finite)}


def bootstrap_summary(predictors, groups, labels, decisions, selected_native, replicates=BOOTSTRAPS,
                      donor_universe=None, group_universe=None):
    samples = {name: {m: [] for m in ("coverage", "error_per_eligible", "error_per_accepted")}
               for name in decisions}
    differences = {mode: {m: [] for m in ("error_per_eligible", "error_per_accepted")}
                   for mode in ("cutoff", "matched_0.95") if f"{mode}:logistic" in decisions}
    for w in bootstrap_weights(predictors, groups, replicates, SEED, donor_universe, group_universe):
        measurements = {name: metric(labels, accept, w) for name, accept in decisions.items()}
        for name, result in measurements.items():
            for m in samples[name]:
                samples[name][m].append(result[m])
        for mode in differences:
            left, right = measurements[f"{mode}:logistic"], measurements[f"{mode}:{selected_native}"]
            for m in differences[mode]:
                differences[mode][m].append(left[m] - right[m] if left[m] is not None and right[m] is not None else None)
    return {"method": "percentile two-way donor x positive-overlap locus-group bootstrap",
            "replicates": replicates, "seed": SEED,
            "conditional_on_frozen_model_rankings_and_decisions": True,
            "intervals": {name: {m: percentile_interval(v) for m, v in ms.items()} for name, ms in samples.items()},
            "paired_logistic_minus_selected_native":
                {mode: {m: percentile_interval(v) for m, v in ms.items()} for mode, ms in differences.items()}}


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, (np.integer, np.bool_)):
        return value.item()
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


def write_json(path, value):
    path.write_text(json.dumps(clean_json(value), indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_csv(path, rows):
    if not rows:
        return
    with gzip.GzipFile(filename="", mode="wb", fileobj=path.open("wb"), mtime=0) as binary:
        import io
        with io.TextIOWrapper(binary, encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)


def read_cosigt(path: Path):
    """Read evaluation-only tables; only byte-identical duplicates may collapse."""
    rows, raw_rows, duplicate_keys = {}, {}, []
    with gzip.open(path, "rt", newline="") as fh:
        header = next(csv.reader([next(fh)], delimiter="\t"))
        for line in fh:
            values = next(csv.reader([line], delimiter="\t"))
            if len(values) != len(header):
                raise ValueError(f"Malformed COSIGT row in {path}")
            row = dict(zip(header, values))
            key = (row["sample"], row["region"])
            if key in rows:
                if raw_rows[key] != line:
                    raise ValueError(f"Non-identical duplicate COSIGT key in {path}: {key}")
                duplicate_keys.append(key)
                continue
            raw_rows[key], rows[key] = line, row
    return rows, {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                  "unique_keys": len(rows), "byte_identical_duplicates_removed": len(duplicate_keys),
                  "duplicate_keys": duplicate_keys}


def describe_cosigt_values(values):
    values = np.asarray(values)
    valid = np.isfinite(values) & (values >= 0)
    return {"n": len(values), "scorable": int(valid.sum()), "unscorable": int((~valid).sum()),
            "continuous": quantiles(values[valid]),
            "discrepancy_counts": {str(t): int((values[valid] > t).sum()) for t in ERROR_THRESHOLDS}}


def cosigt_audit(data_dir: Path, output: Path):
    report = {"scope": "Secondary descriptive reanalysis only; all error/QV fields require truth",
              "confidence_filtering_performed": False,
              "metrics_differ_from_primary": "Published average of two haplotype error rates, not the primary length-weighted diploid div",
              "cohorts": {}}
    paired_rows = []
    for cohort in ("cmrg", "sv"):
        full, full_meta = read_cosigt(data_dir / f"{cohort}_full_30x.tsv.gz")
        all_out, all_meta = read_cosigt(data_dir / f"{cohort}_all_out_30x.tsv.gz")
        common = sorted(set(full) & set(all_out))
        arrays = {
            "full_cosigt": np.asarray([number(full[k]["avg_error_rate_cosigt"]) for k in common]),
            "full_locityper": np.asarray([number(full[k]["avg_error_rate_locityper"]) for k in common]),
            "all_out_cosigt": np.asarray([number(all_out[k]["avg_error_rate_pred"]) for k in common]),
            "all_out_best_available_comparator": np.asarray([number(all_out[k]["avg_error_rate_best"]) for k in common]),
        }
        paired_difference = {}
        for name, first, second in (
                ("cosigt_all_out_minus_full", "all_out_cosigt", "full_cosigt"),
                ("cosigt_all_out_minus_best", "all_out_cosigt", "all_out_best_available_comparator"),
                ("full_cosigt_minus_locityper", "full_cosigt", "full_locityper")):
            a, b = arrays[first], arrays[second]
            valid = np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b >= 0)
            paired_difference[name] = {
                "continuous_difference": quantiles(a[valid] - b[valid]),
                "first_error_second_nonerror": {str(t): int(((a[valid] > t) & (b[valid] <= t)).sum()) for t in ERROR_THRESHOLDS},
                "first_nonerror_second_error": {str(t): int(((a[valid] <= t) & (b[valid] > t)).sum()) for t in ERROR_THRESHOLDS},
            }
        report["cohorts"][cohort] = {
            "full_input": full_meta, "all_out_input": all_meta,
            "matched_keys": len(common), "matched_donors": len({k[0] for k in common}),
            "matched_regions": len({k[1] for k in common}),
            "full_only_keys": sorted(set(full) - set(all_out)),
            "all_out_only_keys": sorted(set(all_out) - set(full)),
            "matched_descriptive": {name: describe_cosigt_values(value) for name, value in arrays.items()},
            "paired_differences": paired_difference,
        }
        for i, (sample, region) in enumerate(common):
            paired_rows.append({"cohort": cohort, "sample": sample, "region": region,
                                **{name: float(value[i]) for name, value in arrays.items()}})
    write_json(output / "cosigt_descriptive.json", report)
    write_csv(output / "cosigt_paired_truth_metrics.csv.gz", paired_rows)
    return report


def run(data_dir: Path, output: Path, groups_path: Path, replicates=BOOTSTRAPS):
    output.mkdir(parents=True, exist_ok=True)
    full_p, full_t, full_meta = read_locityper(data_dir / "full_illumina.csv.gz")
    loo_p, loo_t, loo_meta = read_locityper(data_dir / "loo_illumina.csv.gz")
    if set(full_p) != set(loo_p):
        raise ValueError("Expected the full and LOO tables to contain the identical genotype keys")
    donors = {p.sample for p in loo_p.values()}
    loci = {p.locus for p in loo_p.values()}
    if len(loo_p) != 13040 or len(donors) != 40 or len(loci) != 326:
        raise ValueError("Frozen study requires 13,040 genotype keys, 40 donors, and 326 loci")
    if set(loo_p) != {(s, l) for s in donors for l in loci}:
        raise ValueError("Frozen study requires the complete donor-by-locus matrix")
    with groups_path.open() as fh:
        group_rows = list(csv.DictReader(fh, delimiter="\t"))
    groups = {r["locus"]: r["locus_group"] for r in group_rows}
    loci = {p.locus for p in loo_p.values()}
    if set(groups) != loci or len(groups) != len(group_rows):
        raise ValueError("Locus overlap-group map must cover all loci exactly once")
    if len(set(groups.values())) != 317:
        raise ValueError("Frozen study requires 317 locus overlap groups")
    donor_split = split_ids([p.sample for p in loo_p.values()], "donor")
    group_split = split_ids(list(groups.values()), "locus")
    locus_split = {locus: group_split[group] for locus, group in groups.items()}
    for name, column, expected in (("donor_splits.tsv", "sample", donor_split),
                                  ("locus_splits.tsv", "locus", locus_split)):
        with (groups_path.parent / name).open() as fh:
            frozen_rows = list(csv.DictReader(fh, delimiter="\t"))
        observed = {r[column]: r["split"] for r in frozen_rows}
        if observed != expected or len(observed) != len(frozen_rows):
            raise ValueError(f"Assignments differ from frozen {name}")
        for r in frozen_rows:
            kind, value = ("donor", r["sample"]) if column == "sample" else ("locus", groups[r["locus"]])
            if column == "locus" and r["locus_group"] != groups[r["locus"]]:
                raise ValueError("Frozen locus group differs from provided map")
            if r["assignment_hash"] != stable_hash(kind, value):
                raise ValueError(f"Incorrect assignment hash in frozen {name}")
    keys = sorted(loo_p)
    def category(key):
        return donor_split[key[0]] if donor_split[key[0]] == locus_split[key[1]] else "off_diagonal_unused"
    split_manifest = [{"sample": s, "locus": l, "locus_group": groups[l],
                       "donor_split": donor_split[s], "locus_split": locus_split[l],
                       "analysis_split": category((s, l))} for s, l in keys]
    write_json(output / "split_manifest.json", {"donors": donor_split, "overlap_groups": group_split,
                                               "loci": locus_split, "cases": split_manifest})
    train_keys = [k for k in keys if category(k) == "development" and loo_t[k].scorable]
    cal_keys = [k for k in keys if category(k) == "calibration" and loo_t[k].scorable]
    test_keys = [k for k in keys if category(k) == "test" and loo_t[k].scorable]
    train = [loo_p[k] for k in train_keys]
    cal = [loo_p[k] for k in cal_keys]
    test = [loo_p[k] for k in test_keys]
    y_train = np.asarray([loo_t[k].div > .005 for k in train_keys])
    y_cal = np.asarray([loo_t[k].div > .005 for k in cal_keys])
    model = FrozenLogistic().fit(train, y_train)
    write_json(output / "model.json", model.metadata())
    cal_scores = scores(cal, model)
    cutoffs = {name: calibration_cutoff(risk, cal) for name, risk in cal_scores.items() if name != "no_filter"}
    calibration = {name: metric(y_cal, matched_accept(risk, cal, .95)) for name, risk in cal_scores.items()
                   if name != "no_filter"}
    selected_native = min(NATIVE_METHODS, key=lambda name: (calibration[name]["accepted_wrong"], NATIVE_METHODS.index(name)))
    calibration["no_filter"] = metric(y_cal, np.ones(len(cal), dtype=bool))
    metadata = {"full_input": full_meta, "loo_input": loo_meta,
                "groups_sha256": hashlib.sha256(groups_path.read_bytes()).hexdigest(),
                "analysis_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                            "scipy": scipy.__version__, "scikit_learn": sklearn.__version__},
                "scope": "Retrospective analysis of archived outputs; not rerun alignments or novel algorithm",
                "primary_error": "gt div > 0.005 among finite nonnegative div and finite positive size", "strict_inequality": True,
                "feature_whitelist": ["quality", "weight_dist", "unexpl_reads", "total_reads"],
                "never_features": list(TRUTH_NAMES) + ["sample", "locus", "genotype", "panel_condition"],
                "warnings_unused": True,
                "split_counts": {split: {"all_cases": sum(category(k) == split for k in keys),
                                           "scorable_loo_truth": sum(category(k) == split and loo_t[k].scorable for k in keys),
                                           "scorable_full_truth": sum(category(k) == split and full_t[k].scorable for k in keys),
                                           "unscorable_loo_truth": sum(category(k) == split and not loo_t[k].scorable for k in keys),
                                           "unscorable_full_truth": sum(category(k) == split and not full_t[k].scorable for k in keys)}
                                 for split in ("development", "calibration", "test", "off_diagonal_unused")},
                "selected_native": selected_native, "selection_criterion": "fewest calibration errors retained at ceil(0.95*N); fixed native-order ties",
                "cutoffs": {name: {"risk": value[0], "risk_nonfinite": not np.isfinite(value[0]), "tie_hash": value[1]}
                            for name, value in cutoffs.items()},
                "calibration": calibration,
                "all_cohort_descriptive": {"loo": describe(loo_t), "full": describe(full_t)}}
    all_predictor_rows, all_truth_rows = [], []
    conditions = {}
    for condition, pp, tt in (("loo", loo_p, loo_t), ("full", full_p, full_t)):
        condition_keys = [k for k in test_keys if tt[k].scorable]
        predictors = [pp[k] for k in condition_keys]
        risk_scores = scores(predictors, model)
        decisions = {}
        for name, risk in risk_scores.items():
            if name == "no_filter":
                decisions["cutoff:no_filter"] = np.ones(len(predictors), dtype=bool)
                continue
            decisions[f"cutoff:{name}"] = threshold_accept(risk, predictors, cutoffs[name])
            for retention in RETENTIONS:
                decisions[f"matched_{retention}:{name}"] = matched_accept(risk, predictors, retention)
        sensitivity = {}
        for threshold in ERROR_THRESHOLDS:
            labels = np.asarray([tt[k].div > threshold for k in condition_keys])
            sensitivity[str(threshold)] = {name: metric(labels, accept) for name, accept in decisions.items()}
        labels = np.asarray([tt[k].div > .005 for k in condition_keys])
        conditions[condition] = {
            "n": len(predictors), "truth_missing_excluded_from_loo_test_keys": len(test_keys) - len(condition_keys),
            "descriptive": describe({k: tt[k] for k in condition_keys}),
            "metrics_by_divergence_threshold": sensitivity,
            "bootstrap_primary": bootstrap_summary(predictors, groups, labels, decisions, selected_native, replicates,
                                                     [s for s, split in donor_split.items() if split == "test"],
                                                     [g for g, split in group_split.items() if split == "test"]),
            "calibration_cutoffs_transferred_unchanged": True,
            "matched_retention_is_retrospective_batch_ranking_not_deployment_threshold": True,
        }
        for i, p in enumerate(predictors):
            row = asdict(p)
            row.update(condition=condition, locus_group=groups[p.locus], tie_hash=p.tie,
                       unexpl_fraction=p.unexpl_fraction)
            row.update({"risk_" + name: float(v[i]) for name, v in risk_scores.items()})
            row.update({"accept_" + name: int(v[i]) for name, v in decisions.items()})
            all_predictor_rows.append(row)
            truth_row = asdict(tt[p.key])
            truth_row["condition"] = condition
            all_truth_rows.append(truth_row)
    metadata["test"] = conditions
    metadata["full_control_pairing"] = {"loo_finite_test_cases": len(test_keys),
                                         "paired_finite_full_cases": conditions["full"]["n"]}
    write_json(output / "summary.json", metadata)
    write_csv(output / "heldout_predictions.csv.gz", all_predictor_rows)
    write_csv(output / "heldout_truth.csv.gz", all_truth_rows)
    print(json.dumps({"selected_native": selected_native, "output": str(output),
                      "test_n": len(test), "complete": True}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data/panel_adequacy/locityper"))
    parser.add_argument("--groups", type=Path, default=Path("data/panel_adequacy/locus_groups.tsv"))
    parser.add_argument("--output", type=Path, default=Path("results/panel_adequacy_baseline"))
    parser.add_argument("--cosigt-dir", type=Path, default=Path("data/panel_adequacy/cosigt"))
    args = parser.parse_args()
    run(args.data_dir, args.output, args.groups)
    cosigt_audit(args.cosigt_dir, args.output)
