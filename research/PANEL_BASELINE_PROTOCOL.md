# Reference-panel adequacy: frozen historical baseline study

Protocol written 9 October 2026, before calculating new aggregate outcomes or
fitting filters. This is a pre-specified **retrospective reanalysis**: published
paper summaries were already known. It is not an independently preregistered
new sequencing experiment. The Git commit containing this protocol and input
manifest is the freeze point; its hash will be recorded with the results.

## Question and achievable scope

Do ordinary, deployable confidence diagnostics leave substantial sequence
discrepancies among retained calls when the filter is evaluated on donors and
loci absent from its development and calibration sets?

The primary experiment uses archived Locityper Illumina leave-one-out calls,
with matched full-panel controls. It evaluates **confidence filtering**, not a
new caller or alignment algorithm. It does not yet establish that residual
errors arise from missing structural alleles or are distinguishable using reads.
The eventual 50–1,000 bp indel study in `NEXT_QUESTION.md` requires independently
curated structural truth and read-level observability; those are not supplied
by these summary tables. This bounded study is the first feasibility gate.

## Frozen inputs and provenance

Primary source: Timofey Prodanov, Locityper benchmark archive, pinned Zenodo
record [15879519](https://doi.org/10.5281/zenodo.15879519), CC-BY-4.0.
Archive SHA256:
`fd0e8ab023d8173c1cacf7627f07f7cacc62bcea5a2ab5f171d096bcc090312f`.
Use `QV/loo_illumina.csv.gz` and `QV/full_illumina.csv.gz`, preserved byte-for-byte
under `data/panel_adequacy/locityper/`. They are tab-separated, with comment
lines. Select **only `query_type == "gt"`**: the two haplotype rows per case are
not extra observations. Inventory found 13,040 unique genotype cases per panel,
40 donors × 326 loci, with identical keys. Verify this on every run.

Version uncertainty is preserved: the archive README says v0.17.3; embedded
commands name v0.17.0 input summaries. These are the archived results, not a
new run of either version or a claim about current Locityper performance.
The upstream evaluation script omits loci without comparison files: the
denominator is **published evaluable calls**, not all attempted calls.

Secondary source: the COSIGT paper repository at commit
`3d84bf8490bbff5ec33cf72706e713d1a41730ed`, with its licence and exact source
paths/checksums in `data/panel_adequacy/SOURCE_MANIFEST.json`. Four tables cover
full-panel and cross-panel evaluations for CMRG and large-SV cohorts. Their
truth-derived QV values are not deployable confidence scores. These tables
support a descriptive discrepancy audit only, never confidence-filter training.

## Independence and fixed splits

Use the original COSIGT-bundled `cmrgs.hprcy1.bed`, adapted from the Locityper
targets, as coordinate metadata. It matches all 326 locus names. Form connected
components of positive-length overlap on the same chromosome using half-open
intervals; touching intervals do not overlap. Inventory gives 317 components:
308 singletons and nine pairs. A component ID is its lexically sorted locus
names joined by commas. Preserve the complete map in `locus_groups.tsv`.

Independently order donor IDs and locus component IDs by hexadecimal SHA256 of:

* `panel-native-v1|donor|` followed by the donor ID;
* `panel-native-v1|locus|` followed by the component ID.

For each ordering allocate the first `floor(n/2)` identities to development,
the next `floor(n/4)` to calibration, and the remainder to test. This gives
20/10/10 donors and 158/79/80 locus components. The outcome-blind split manifest
gives 161/84/81 loci and 3,220/840/810 diagonal cases respectively. All members
of a component inherit its split.

Fit only on development-donor × development-locus cases. Select scalar filters
and retention thresholds only on calibration × calibration. Evaluate the
confirmatory comparison once on test × test. Off-diagonal cases are unused
for fitting/selection and appear only in explicitly descriptive cohort totals.
Save all identity assignments before outcome analysis.

Donor relatives and non-overlapping homologous loci are not grouped: suitable
curated metadata was not bundled. Consequently this is donor-ID and overlapping
interval separation, not a guarantee of ancestry, family or gene-family
independence. The original caller/panel construction predates these splits;
only our post-hoc filter development obeys them.

## Outcomes and exclusions

Primary error proxy: archived genotype `div > 0.005` (more than 0.5% sequence
discrepancy). Predeclared sensitivities: `div > 0.001` and `div > 0.01`. Report
continuous divergence as context. `div` is the diploid comparison's sum of edit
counts divided by sum of alignment lengths, after the upstream pairing rule.
It is not a structural-genotype label or a probability that a call is wrong.

Require finite, nonnegative `div` and positive comparison `size`; report
unscorable counts by panel and split, and exclude them from outcome rates
without silently treating them as correct. Do not exclude legitimate infinite
`qv` values. Preserve finite native features independently from truth validity.
If the development outcome has fewer than two classes, do not fit logistic
regression; report the combined baseline as unavailable, rather than change the
outcome or tune the split. Unexpected schema/key failures stop the run.

LOO `avail_div` is the **published closest-available comparator**. Its upstream
procedure selects each truth haplotype's closest available sequence independently
and aggregates their distances; with varying alignment lengths this need not
minimize the diploid ratio globally. Do not call it a certified error floor.
Use it only to describe agreement and discrepancy, not as a predictor or to
declare an error biologically avoidable.

## Frozen predictors and ordinary baselines

The input feature whitelist is `quality`, `weight_dist`, `unexpl_reads`,
`total_reads`, and `warnings`. Warnings are constant `*` in the primary archive
and therefore provide no discrimination. No label fields (`edit`, `size`,
`div`, `qv`, any `avail_*`), haplotype IDs, genotype IDs or donor/locus identities
enter a statistical predictor. In particular, truth alignment `size` cannot be
used to normalize coverage. Identities are used only for splits, ties and audit.

Scalar risk scores, lower being safer, in fixed tie-preference order:

1. Negative genotype quality.
2. Weighted distance.
3. Unexplained-read fraction (`unexpl_reads / total_reads`).
4. Unexplained-read count.
5. Negative total recruited-read count.

These are general native-diagnostic baselines. Do not describe the count or
coverage filters as an author-recommended general filtering policy.
Missing scalar risks sort last (positive infinity); no outcome-dependent
imputation is permitted. Require nonnegative native numeric inputs before
logarithmic transforms; unexpected negative values stop the run.

The ordinary combined baseline is fixed L2 logistic regression (`C=1`, no class
weight, maximum 2,000 iterations, random seed 20261009) on five features:
`log1p(quality)`, `log1p(weight_dist)`, `log1p(unexpl_reads)`,
`log1p(total_reads)`, and unexplained-read fraction. Zero-denominator fractions
are missing. Use development-only median imputation with missing indicators
and development-only standardization. No feature search, hyperparameter tuning
or outcome-dependent transformations. Record convergence and coefficients.

Select the strongest **scalar** baseline by its number of primary-error calls
retained among the safest `ceil(0.95*N)` calibration cases. Break equal error
counts by the fixed list above. Evaluate the logistic baseline separately;
do not choose a winner between it and scalar baselines using test outcomes.
Include hash-only selection as a descriptive outcome-independent control.

## Retention, ties and transferred thresholds

Break score ties by SHA256 of `panel-native-v1|tie|<sample>|<locus>` (lexical
hexadecimal order). This cannot use outcomes. For each score freeze the tuple
`(risk, tie_hash)` of the `ceil(0.95*N)`th calibration case as an inclusive
threshold. Transfer that tuple unchanged to test and report achieved retention,
which need not equal 95%.

Separately report retrospective matched-retention rankings at 80%, 90%, 95%,
99% and 100%, using `ceil(retention*N)` and the same ties. These compare
discrimination at the same test call budget but are not evidence that a fixed
deployed threshold achieves that retention. Report exact numerator/denominator
counts and actual fractions. Reapply the frozen LOO model and thresholds to
the identical full-panel test keys without refitting. Matched-retention full
controls are ranked separately and identified as such.

Report retained error count, error fraction among accepted calls, retained
error count per eligible case, and fraction of baseline errors removed. For
sensitivities, reuse the primary fitted model, scalar selection and decisions;
do not retrain for a more favourable endpoint.

## Uncertainty and interpretation

Use 4,999 paired two-way bootstrap replicates, seed 20261009, resampling the
ten held-out donor IDs and held-out locus components independently. Case weights
are products of donor and component multiplicities; use the same weights for
all methods and both panels. Report percentile 95% intervals for retained error
fractions and paired logistic-minus-selected-scalar error differences per
eligible case. Decisions and model parameters remain fixed: these intervals
are conditional on the fitted filters and observed cohorts, not uncertainty
over retraining or undisclosed family/homology relationships. No row-independent
confidence intervals. Few errors or wide intervals are inconclusive.

This study establishes whether a residual sequence-discrepancy problem exists
after ordinary filters; it cannot establish a new algorithm's advantage.
Do not apply the later 20% structural-error improvement gate to these different
labels. A positive headroom result warrants a small read-level audit; it does
not warrant immediate development of a novel aligner. If errors are scarce,
explained by ordinary filters, or unidentifiable from reads, stop the proposed
algorithm. Genotype–phenotype benefit remains untested, and selective missingness
could introduce association bias.

## Secondary descriptive and engineering work

For COSIGT preserve both cohorts separately, exactly deduplicate byte-identical
rows before pairing (inventory found 24 such rows in the all-out SV table),
report unmatched keys, and summarize published prediction and best-available
sequence comparisons. These cohorts, panel exclusions and metrics differ from
the primary study; do not pool them or rank tools using unmatched samples.
No claims about detecting 50–1,000 bp missing indels follow from the selected
large-SV cohort or its extended loci.

Current Locityper binaries and official two-locus demo inputs may be prepared
as an engineering follow-up. Version checks are not a genotype benchmark.
`PANEL_RUNTIME_FEASIBILITY.md` records exactly what was executed; no current
performance claim is made without a separately specified and completed run.

## Reproduction and audit

Freeze source checksums, licence notices, split assignments and this protocol
on GitHub before the first analysis run. Tests may use synthetic data before
the freeze. Preserve code, environment versions, command, source/commit hashes,
per-case held-out decisions, summaries and independent audit notes. Subsequent
corrections must be identified as corrections, never silently described as
pre-specified. The study report must distinguish observations, inferences and
questions that these data cannot answer.
