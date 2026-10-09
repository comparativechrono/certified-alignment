# Independent audit of the historical baseline study

Status: **passed**. This audit independently reconstructed the reported results
from the archived source tables, without importing the analysis program,
refitting a model, changing splits, or selecting a different endpoint.

The protocol and analysis were frozen at commit
`d6d3b76ff823c85dafed3704667ad1713df17750`. The audited working copies of the
protocol, analysis, donor assignments and locus assignments exactly match that
commit. All copied source checksums and the derived overlap-map checksum match
the source manifest. The study remains a retrospective analysis of archived
calls with unresolved v0.17.0/v0.17.3 provenance; it is not a current-caller
benchmark or a structural-allele validation experiment.

## Verification performed

Before outcomes were calculated, a separate pairwise-interval implementation
checked all 52,975 locus pairs and independently reconstructed 317 overlap
components. Donor and component hashes reproduce the frozen assignments:

| Split | Donors | Components | Loci | Donor–locus cases |
| --- | ---: | ---: | ---: | ---: |
| Development | 20 | 158 | 161 | 3,220 |
| Calibration | 10 | 79 | 84 | 840 |
| Test | 10 | 80 | 81 | 810 |

The implementation review confirmed that only development labels enter the
ordinary logistic regression; only calibration labels select the scalar
baseline. Native feature construction excludes all truth comparison fields,
including alignment length. The audit independently recomputed development-only
medians, means and scales, then manually reconstructed logistic probabilities
from the saved coefficients and native features. It independently verified:

- Calibration selection of **unexplained-read fraction**, including fixed tie
  preference, tuple cutoffs, and every calibration metric.
- Every held-out score and acceptance decision, all three discrepancy
  thresholds and all declared retention budgets, for both panels.
- Main selected-scalar and logistic bootstrap intervals and paired differences,
  using an independent outer-product cluster-weight calculation. These remain
  conditional intervals for the fixed fitted models and decisions.
- COSIGT key matching, duplicate counts, discrepancy counts and metric units.

The executable audit completed 402,514 assertions. Most are per-record numeric
checks; this count is not a measure of statistical power. The source analysis
also passed 12 synthetic tests covering splits, missingness, calibration, ties,
feature isolation, clustered resampling and duplicate handling.

No outcome-run correction or rerun was needed. The analysis was run once after
the freeze; the separate audit was run again after adding the explicitly
post-hoc descriptions below.

## Confirmed primary findings and denominators

Four of the 810 held-out cases are unscorable, leaving **806 evaluable cases**
in each panel. Across the entire archived matrix, 171 of 13,040 cases per panel
are unscorable. Infinite QV for exact agreement was not excluded.

The primary label is `div > 0.005`, a sequence-discrepancy proxy. In the LOO
test set, 83 of 806 cases exceed this threshold. At retrospective 95% retention,
both compared filters retain `ceil(0.95 × 806) = 766` cases:

| Frozen filter | Retained cases | Retained discrepancies |
| --- | ---: | ---: |
| Calibration-selected native scalar | 766 | 63 |
| Ordinary combined logistic baseline | 766 | 65 |

Transferring calibration thresholds unchanged produces different retentions:
the selected native scalar retains 771 cases with 64 discrepancies, while the
logistic baseline retains 767 with 65. These transferred thresholds must not be
reported as exact 95% retention on the test cohort.

## Additional post-hoc descriptions

The following diagnostics were requested after results existed. They do not
change any filter, endpoint, split, exclusion or conclusion test.

The 63 discrepancies retained by the native scalar occur across all ten test
donors and 29 loci. MUC1 contributes ten and TPO six; MUC19, SLC22A1, SMOC2 and
CTDP1 contribute four each. Thus the six leading loci account for 32 of 63
discrepancies; this is appreciable concentration without a single-locus
explanation. These loci were not used to construct a more favourable test set.

For those same 63 cases:

| Published closest-available `avail_div` | Full-panel `div ≤0.005` | Full-panel `div >0.005` | Total |
| --- | ---: | ---: | ---: |
| `≤0.005` | 27 | 10 | 37 |
| `>0.005` | 22 | 4 | 26 |
| Total | 49 | 14 | 63 |

No comparator value is missing in this cross-tabulation. Neither its rows nor
columns establish structural allele absence, reliable truth, read-level
observability, or biological avoidability. The official
[v0.17.3 evaluation code](https://github.com/tprodanov/locityper/blob/146cfd42e9179bafd63d1f41db7fc449940989de/extra/gt_dist.py)
selects nearest sequences per haplotype before aggregating distances. Therefore
`avail_div` need not be the global minimum diploid distance ratio.

There is also a simple label-aware bound on filtering headroom. Retaining 766
of 806 cases permits only 40 rejections. Even a filter given perfect error
labels must retain at least `83 − 40 = 43` discrepancies, or 5.61% of accepted
calls. Relative to the native filter's 63, at most 20 further discrepancies
could be removed at this call budget, a 31.75% relative reduction. This is an
oracle ranking bound, not an achievable method or evidence that reads identify
the right cases. Residual errors alone do not imply arbitrarily large useful
filtering headroom.

## Secondary COSIGT checks

The `avg_error_rate_*` fields are fractions, not percentages. For each finite
entry, the audit checked the mean of the two component error rates and checked
each component against `10^(-QV/10)`; the largest relative inversion difference
was below `2.3 × 10^-16`. All are truth-sequence comparisons, not native
confidence scores. Their equal-weight mean across haplotypes differs from the
primary study's length-weighted diploid divergence.

The matched descriptive cohorts contain 11,674 CMRG cases and 14,838 SV cases.
The 24 identical SV duplicate records were removed once. All matched COSIGT
prediction and best-comparator values are scorable; companion full-panel
Locityper values have five missing CMRG evaluations and 22 missing SV
evaluations. The verified counts above 0.005 are:

| Metric | CMRG discrepancies / scorable cases | SV discrepancies / scorable cases |
| --- | ---: | ---: |
| Full-panel COSIGT | 326 / 11,674 | 585 / 14,838 |
| Full-panel Locityper companion results | 138 / 11,669 | 277 / 14,816 |
| Cross-panel COSIGT | 1,674 / 11,674 | 2,326 / 14,838 |
| Published best available comparator | 629 / 11,674 | 898 / 14,838 |

These are full-cohort descriptive counts, separate from the 806-case held-out
filter experiment. They do not support a matched confidence-filter comparison
between COSIGT and Locityper, nor a claim about 50–1,000 bp indel genotypes.

## Reproduction

From the repository root, after the frozen analysis has produced its outputs:

```sh
python studies/panel_adequacy/audit_baseline.py
```

The audit writes `results/panel_adequacy_baseline/independent_audit.json`, which
contains check counts and the complete descriptive concentration tables. Its
file reads are limited to the source snapshot, frozen git objects, saved model
and saved outputs. It neither downloads new evidence nor trains a predictor.
