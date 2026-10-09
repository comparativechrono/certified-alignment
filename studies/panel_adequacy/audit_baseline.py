#!/usr/bin/env python3
"""Independent evidence audit; does not import analyze.py or refit any model.

Run after the frozen analysis. Reconstructs input identities, transforms, scores,
decisions, metrics and selected bootstrap intervals from the archived bytes.
Additional error concentration summaries are explicitly post-hoc descriptions.
"""
import csv
import gzip
import hashlib
import itertools
import json
import math
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / 'data/panel_adequacy'
OUT = ROOT / 'results/panel_adequacy_baseline'
NATIVE = ['quality', 'weight_dist', 'unexpl_fraction', 'unexpl_reads', 'total_reads']
THRESHOLDS = [.005, .001, .01]
RETENTIONS = [.8, .9, .95, .99, 1.]
CHECKS = Counter()


def check(condition, name):
    if not condition:
        raise AssertionError(name)
    CHECKS[name] += 1


def near(a, b, name, atol=1e-12):
    check(bool(np.allclose(a, b, rtol=1e-10, atol=atol, equal_nan=True)), name)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def h(kind, key):
    return hashlib.sha256(f'panel-native-v1|{kind}|{key}'.encode()).hexdigest()


def tab(path, delimiter='\t'):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader((line for line in f if not line.startswith('#')), delimiter=delimiter))


def num(row, field):
    try:
        return float(row[field])
    except (KeyError, ValueError, TypeError):
        return math.nan


def valid(row):
    return math.isfinite(num(row, 'div')) and num(row, 'div') >= 0 and math.isfinite(num(row, 'size')) and num(row, 'size') > 0


def identity(row):
    return row['sample'], row['locus']


def load_primary(name):
    rows = [r for r in tab(DATA / 'locityper' / f'{name}_illumina.csv.gz') if r['query_type'] == 'gt']
    index = {identity(r): r for r in rows}
    check(len(index) == len(rows) == 13040, 'primary_unique_keys')
    return index


def independently_split(ids, kind):
    ids = sorted(set(ids), key=lambda value: h(kind, value))
    a, b = len(ids) // 2, len(ids) // 4
    return {v: 'development' if i < a else 'calibration' if i < a+b else 'test' for i, v in enumerate(ids)}


def verify_identities(primary):
    bed = [(n, c, int(a), int(b)) for c, a, b, n in (line.split()[:4] for line in (DATA / 'cosigt/cmrgs_original.bed').read_text().splitlines())]
    parent = {n: n for n, *_ in bed}
    def root(n):
        while parent[n] != n:
            n = parent[n]
        return n
    for (n, c, a, b), (m, d, x, y) in itertools.combinations(bed, 2):
        if c == d and max(a, x) < min(b, y):
            parent[root(m)] = root(n)
    components = defaultdict(list)
    for n, *_ in bed:
        components[root(n)].append(n)
    groups = {n: ','.join(sorted(members)) for members in components.values() for n in members}
    check(len(components) == 317, 'independent_overlap_components')
    observed = tab(DATA / 'locus_groups.tsv')
    check({r['locus']: r['locus_group'] for r in observed} == groups, 'overlap_memberships')
    donor_ids = {k[0] for k in primary}
    check(set(primary) == set(itertools.product(donor_ids, groups)), 'complete_donor_locus_matrix')
    donors = independently_split(donor_ids, 'donor')
    group_splits = independently_split(groups.values(), 'locus')
    locus_splits = {l: group_splits[g] for l, g in groups.items()}
    check({r['sample']: r['split'] for r in tab(DATA / 'donor_splits.tsv')} == donors, 'donor_split')
    check({r['locus']: r['split'] for r in tab(DATA / 'locus_splits.tsv')} == locus_splits, 'locus_split')
    split_keys = {s: sorted(k for k in primary if donors[k[0]] == locus_splits[k[1]] == s) for s in ['development', 'calibration', 'test']}
    return groups, donors, group_splits, split_keys


def vector(row):
    raw = [num(row, x) for x in ['quality', 'weight_dist', 'unexpl_reads', 'total_reads']]
    x = [math.log1p(v) if math.isfinite(v) and v >= 0 else math.nan for v in raw]
    x.append(raw[2] / raw[3] if math.isfinite(raw[2]) and raw[3] > 0 else math.nan)
    return np.array(x)


def get_scores(row, model):
    q, w, u, n = [num(row, f) for f in ['quality', 'weight_dist', 'unexpl_reads', 'total_reads']]
    values = [-q, w, u / n if n > 0 else math.inf, u, -n]
    result = {name: x if math.isfinite(x) else math.inf for name, x in zip(NATIVE, values)}
    result.update(hash_control=0., no_filter=0.)
    if model['available']:
        x = vector(row)
        missing = ~np.isfinite(x)
        z = np.concatenate([np.where(missing, model['imputation_medians'], x), missing.astype(float)])
        z = (z - model['scaler_mean']) / model['scaler_scale']
        logit = float(z @ np.array(model['coefficients_standardized']) + model['intercept'])
        result['logistic'] = 1 / (1 + math.exp(-logit))
    return result


def tie(key):
    return h('tie', '|'.join(key))


def ordered_keys(keys, risks, method):
    return sorted(keys, key=lambda key: (risks[key][method], tie(key)))


def independent_metrics(rows, keys, accepted, threshold):
    n = len(keys)
    wrong = sum(num(rows[k], 'div') > threshold for k in accepted)
    total = sum(num(rows[k], 'div') > threshold for k in keys)
    return dict(eligible=float(n), accepted=float(len(accepted)), accepted_wrong=float(wrong),
                coverage=len(accepted)/n, error_per_eligible=wrong/n,
                error_per_accepted=wrong/len(accepted) if accepted else None,
                total_wrong=float(total), fraction_errors_removed=1-wrong/total if total else None)


def verify_metrics(actual, expected):
    for field, value in expected.items():
        if value is None:
            check(actual[field] is None, 'metric_null')
        else:
            near(actual[field], value, 'primary_metric')


def main():
    summary = json.loads((OUT / 'summary.json').read_text())
    model = json.loads((OUT / 'model.json').read_text())
    provenance = json.loads((OUT / 'run_provenance.json').read_text())
    source_manifest = json.loads((DATA / 'SOURCE_MANIFEST.json').read_text())
    for entry in source_manifest['files'] + source_manifest['derived_files']:
        check(sha(ROOT / entry['path']) == entry['sha256'], 'source_sha256')
    for path in ['research/PANEL_BASELINE_PROTOCOL.md', 'studies/panel_adequacy/analyze.py', 'data/panel_adequacy/donor_splits.tsv', 'data/panel_adequacy/locus_splits.tsv']:
        frozen = subprocess.check_output(['git', 'show', provenance['protocol_freeze_commit'] + ':' + path], cwd=ROOT)
        check(frozen == (ROOT / path).read_bytes(), 'freeze_bytes_unchanged')
    check(sha(ROOT / 'studies/panel_adequacy/analyze.py') == provenance['analysis_sha256'] == summary['analysis_sha256'], 'analysis_code_hash')
    primary = {name: load_primary(name) for name in ['loo', 'full']}
    check(set(primary['loo']) == set(primary['full']), 'paired_primary_keys')
    groups, donors, group_splits, split_keys = verify_identities(primary['loo'])
    train = [primary['loo'][k] for k in split_keys['development'] if valid(primary['loo'][k])]
    if model['available']:
        x = np.array([vector(r) for r in train])
        medians = np.array([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else 0. for v in x.T])
        missing = ~np.isfinite(x)
        z = np.concatenate([np.where(missing, medians, x), missing.astype(float)], axis=1)
        scale = z.std(axis=0); scale[scale == 0] = 1
        near(medians, model['imputation_medians'], 'development_only_imputation')
        near(z.mean(axis=0), model['scaler_mean'], 'development_only_mean')
        near(scale, model['scaler_scale'], 'development_only_scaling')
    cal = [k for k in split_keys['calibration'] if valid(primary['loo'][k])]
    risk_cal = {k: get_scores(primary['loo'][k], model) for k in cal}
    methods = list(risk_cal[cal[0]])
    cal_keep = {m: set(ordered_keys(cal, risk_cal, m)[:math.ceil(.95*len(cal))]) for m in methods if m != 'no_filter'}
    cal_errors = {m: sum(num(primary['loo'][k], 'div') > .005 for k in cal_keep[m]) for m in NATIVE}
    selected = min(NATIVE, key=lambda m: (cal_errors[m], NATIVE.index(m)))
    check(selected == summary['selected_native'], 'calibration_scalar_winner')
    cutoffs = {}
    for m in methods:
        if m == 'no_filter':
            continue
        last = ordered_keys(cal, risk_cal, m)[math.ceil(.95*len(cal))-1]
        cutoffs[m] = risk_cal[last][m], tie(last)
        near(cutoffs[m][0], summary['cutoffs'][m]['risk'], 'calibration_risk_cutoff')
        check(cutoffs[m][1] == summary['cutoffs'][m]['tie_hash'], 'calibration_tie_cutoff')
        verify_metrics(summary['calibration'][m], independent_metrics(primary['loo'], cal, cal_keep[m], .005))
    saved_predictions = {(r['condition'], r['sample'], r['locus']): r for r in tab(OUT / 'heldout_predictions.csv.gz', ',')}
    report = {'status': 'pass', 'freeze_commit': provenance['protocol_freeze_commit'], 'calibration_selected_scalar': selected,
              'scope': 'Independent recomputation; no model refitting or split changes; no import of analyze.py',
              'primary': {}, 'cosigt': {}}
    test_loo = [k for k in split_keys['test'] if valid(primary['loo'][k])]
    donor_order = sorted(d for d, s in donors.items() if s == 'test')
    group_order = sorted(g for g, s in group_splits.items() if s == 'test')
    for condition, rows in primary.items():
        keys = [k for k in test_loo if valid(rows[k])]
        risks = {k: get_scores(rows[k], model) for k in keys}
        decisions = {'cutoff:no_filter': set(keys)}
        for m in methods:
            if m == 'no_filter':
                continue
            decisions['cutoff:'+m] = {k for k in keys if (risks[k][m], tie(k)) <= cutoffs[m]}
            for retention in RETENTIONS:
                decisions[f'matched_{retention}:{m}'] = set(ordered_keys(keys, risks, m)[:math.ceil(retention*len(keys))])
        for k in keys:
            saved = saved_predictions[(condition, *k)]
            for m, value in risks[k].items():
                near(float(saved['risk_'+m]), value, 'manual_native_and_logistic_scores')
            for name, keep in decisions.items():
                check(int(saved['accept_'+name]) == int(k in keep), 'heldout_decision')
        for threshold in THRESHOLDS:
            for name, keep in decisions.items():
                verify_metrics(summary['test'][condition]['metrics_by_divergence_threshold'][str(threshold)][name], independent_metrics(rows, keys, keep, threshold))
        # Independently calculate main paired bootstrap statistics using summed
        # outer-product cluster count matrices rather than production metrics.
        names = [f'{mode}:{m}' for mode in ['cutoff', 'matched_0.95'] for m in [selected, 'logistic'] if m in methods]
        di = np.array([donor_order.index(k[0]) for k in keys])
        gi = np.array([group_order.index(groups[k[1]]) for k in keys])
        labels = np.array([num(rows[k], 'div') > .005 for k in keys])
        masks = np.array([[k in decisions[name] for k in keys] for name in names])
        intervals = defaultdict(list); differences = defaultdict(list)
        rng = np.random.default_rng(20261009)
        for _ in range(4999):
            d = np.bincount(rng.integers(10, size=10), minlength=10)
            g = np.bincount(rng.integers(80, size=80), minlength=80)
            weights = np.outer(d, g)[di, gi]
            retained = masks @ weights
            wrong = (masks & labels) @ weights
            den = weights.sum()
            per_eligible = wrong / den
            per_accepted = np.divide(wrong, retained, out=np.full(len(wrong), np.nan), where=retained != 0)
            for i, name in enumerate(names):
                intervals[name, 'error_per_eligible'].append(per_eligible[i])
                intervals[name, 'error_per_accepted'].append(per_accepted[i])
            if 'logistic' in methods:
                for mode in ['cutoff', 'matched_0.95']:
                    a, b = names.index(mode+':logistic'), names.index(mode+':'+selected)
                    differences[mode, 'error_per_eligible'].append(per_eligible[a]-per_eligible[b])
                    differences[mode, 'error_per_accepted'].append(per_accepted[a]-per_accepted[b])
        bootstrap = summary['test'][condition]['bootstrap_primary']
        for (name, metric), values in intervals.items():
            values = np.asarray(values); values = values[np.isfinite(values)]
            obs = bootstrap['intervals'][name][metric]
            near(np.quantile(values, [.025, .975]), [obs['low'], obs['high']], 'main_bootstrap_intervals')
        for (mode, metric), values in differences.items():
            values = np.asarray(values); values = values[np.isfinite(values)]
            obs = bootstrap['paired_logistic_minus_selected_native'][mode][metric]
            near(np.quantile(values, [.025, .975]), [obs['low'], obs['high']], 'paired_bootstrap_intervals')
        # These concentrations were requested after the freeze, are descriptive,
        # and are never used to remove a locus or choose a better test subset.
        concentrations = {}
        for name in ['cutoff:no_filter', f'matched_0.95:{selected}', 'matched_0.95:logistic']:
            if name not in decisions:
                continue
            bad = [k for k in keys if k in decisions[name] and num(rows[k], 'div') > .005]
            total_errors = sum(num(rows[k], 'div') > .005 for k in keys)
            rejected = len(keys) - len(decisions[name])
            oracle_errors = max(0, total_errors - rejected)
            concentrations[name] = {
                'errors': len(bad), 'donors_with_errors': len({k[0] for k in bad}), 'loci_with_errors': len({k[1] for k in bad}),
                'groups_with_errors': len({groups[k[1]] for k in bad}),
                'errors_by_donor': dict(Counter(k[0] for k in bad).most_common()),
                'errors_by_locus': dict(Counter(k[1] for k in bad).most_common()),
                'closest_available_above_005': sum(num(rows[k], 'avail_div') > .005 for k in bad) if condition == 'loo' else None,
                'closest_available_at_or_below_005': sum(num(rows[k], 'avail_div') <= .005 for k in bad) if condition == 'loo' else None,
                'closest_available_unscorable': sum(not math.isfinite(num(rows[k], 'avail_div')) for k in bad) if condition == 'loo' else None,
                'label_aware_rejection_budget_bound': {
                    'interpretation': 'Post-hoc perfect-error-label ranking bound; not a deployable filter or evidence of observability',
                    'eligible': len(keys), 'accepted': len(decisions[name]), 'rejection_budget': rejected,
                    'initial_errors': total_errors, 'minimum_retained_errors': oracle_errors,
                    'minimum_error_fraction_among_accepted': oracle_errors / len(decisions[name]) if decisions[name] else None,
                    'maximum_further_errors_removable': len(bad) - oracle_errors,
                    'maximum_relative_error_reduction': (len(bad)-oracle_errors)/len(bad) if bad else None,
                },
            }
            if condition == 'loo':
                cross = Counter()
                for key in bad:
                    available = num(rows[key], 'avail_div')
                    a = 'above' if math.isfinite(available) and available > .005 else 'at_or_below' if math.isfinite(available) else 'unscorable'
                    full = primary['full'][key]
                    b = 'above' if valid(full) and num(full, 'div') > .005 else 'at_or_below' if valid(full) else 'unscorable'
                    cross[f'available_{a}|full_control_{b}'] += 1
                concentrations[name]['closest_available_by_full_control_sequence_discrepancy'] = dict(cross)
        report['primary'][condition] = {'heldout_cases': len(keys), 'unscorable_all_source_cases': sum(not valid(r) for r in rows.values()),
                                        'posthoc_descriptive_error_concentration': concentrations}
    secondary = json.loads((OUT / 'cosigt_descriptive.json').read_text())
    for cohort in ['cmrg', 'sv']:
        tables = {}; duplicate_count = {}
        for panel in ['full', 'all_out']:
            path = DATA / 'cosigt' / f'{cohort}_{panel}_30x.tsv.gz'
            records = tab(path); indexed = {}; duplicates = 0
            for row in records:
                key = row['sample'], row['region']
                if key in indexed:
                    check(indexed[key] == row, 'cosigt_duplicate_identical_fields'); duplicates += 1
                indexed[key] = row
            tables[panel] = indexed; duplicate_count[panel] = duplicates
        common = sorted(set(tables['full']) & set(tables['all_out']))
        csum = secondary['cohorts'][cohort]
        check(len(common) == csum['matched_keys'], 'cosigt_matched_keys')
        check(duplicate_count['all_out'] == csum['all_out_input']['byte_identical_duplicates_removed'], 'cosigt_duplicate_count')
        counts = {}; max_qv_relative_error = 0.
        for panel, suffix, name in [('full','cosigt','full_cosigt'), ('full','locityper','full_locityper'), ('all_out','pred','all_out_cosigt'), ('all_out','best','all_out_best_available_comparator')]:
            values = []; qv_checked = 0
            for key in common:
                row = tables[panel][key]
                value = num(row, 'avg_error_rate_'+suffix)
                e = [num(row, f'error_rate_{i}_'+suffix) for i in [1,2]]
                if all(math.isfinite(x) for x in e+[value]):
                    near(value, sum(e)/2, 'cosigt_mean_haplotype_error', atol=1e-8)
                for i in [1,2]:
                    q = num(row, f'QV_{i}_'+suffix)
                    err = num(row, f'error_rate_{i}_'+suffix)
                    if math.isfinite(q) and math.isfinite(err):
                        reconstructed = 10**(-q/10)
                        near(err, reconstructed, 'cosigt_error_is_fraction_from_qv', atol=1e-7)
                        if reconstructed:
                            max_qv_relative_error = max(max_qv_relative_error, abs(err-reconstructed)/reconstructed)
                        qv_checked += 1
                if math.isfinite(value) and value >= 0:
                    values.append(value)
            stats = csum['matched_descriptive'][name]
            check(stats['scorable'] == len(values), 'cosigt_scorable_count')
            actual = {str(t): sum(v > t for v in values) for t in THRESHOLDS}
            check(actual == stats['discrepancy_counts'], 'cosigt_discrepancy_counts')
            counts[name] = {'scorable': len(values), 'unscorable': len(common)-len(values), 'threshold_counts': actual, 'qv_fraction_checks': qv_checked}
        report['cosigt'][cohort] = {'matched_cases': len(common), 'duplicates': duplicate_count, 'unit': 'fractions, not percentages',
                                    'metrics': counts, 'max_qv_inversion_relative_error': max_qv_relative_error}
    report['checks'] = dict(CHECKS)
    (OUT / 'independent_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'status':'pass', 'check_assertions':sum(CHECKS.values()), 'selected_scalar':selected,
                      'primary_test_cases':{c:r['heldout_cases'] for c,r in report['primary'].items()}}))


if __name__ == '__main__':
    main()
