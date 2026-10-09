# Fresh tandem-repeat validation: frozen design

Fixed 9 October 2026 before generating any new alignment result or timing.
This is a validation experiment using the unchanged v4 engine, not a fifth
alignment algorithm. V4's retrospective nine-case tandem-repeat result motivates
the test but contributes no observations to its endpoints.

## Question and fixed candidate

Does **linear/sparse** retain a useful total-runtime advantage on fresh tandem
repeat pairs, relative to **both linear/length and linear/qgram**? The dense
linear/ordered bound is a secondary comparison. The candidate was chosen from
the previous experiment and will not be reselected on these data.

The primary endpoint is the ratio of sums of per-case median internal total
runtime on the 60 fresh direct-replication pairs. Confirmation requires sparse
to be at least 10% faster than each baseline (both ratios <=0.90), both paired
cluster-bootstrap 95% upper limits below 1.00, and all exact-score validations
to pass. Seven timed repetitions and one warmup are fixed before timing. This
is a practical experiment gate, not a claim of family-wise statistical coverage.
No parameter or generator tuning after observing outcomes is allowed.

## Fresh data: 380 pairs

**Direct replication, 60 pairs.** Reuse `run_bench.generate_case` unchanged for
scenario `tandem_repeats`, source lengths 256, 1,024 and 4,096 and new seeds
501..520. This preserves the original roughly 75% repeat tract, shared random
flanks, randomly drawn motif length 2..6, whole-copy expansion at the tract
entrance, roughly 2% length change and 0.4% target substitution process. Motifs
in this replication remain potentially nonprimitive, as in the original code.
All 20 seeds contribute all three lengths. None of these pairs was used in v4.

**Controlled transfer, 160 repeat pairs.** There are 20 independently generated
loci: source lengths {1,024,4,096} x primitive motif periods {2,3,4,5,6} x seeds
{601,607}. The central tract contains `floor((n//2)/period)` whole motif copies;
the remaining bases form random flanks, so source length is exactly n. A motif
must not be a repetition of a shorter motif. Each locus has eight variants:

- Expansion or contraction at a whole-copy boundary near the tract midpoint.
- Gap size `period * max(1, floor(fraction*n/period))`, fraction 0.02 or 0.08.
- Perfect regime: no interruptions or target substitutions. Imperfect regime:
  1% source-tract substitution probability, inherited through the edit, and
  0.4% target substitution probability globally.

Expansion inserts exact motif copies; contraction removes the actual source
substring, including any interruptions. The imperfect regime changes reference
purity and pair noise together; its effect cannot be attributed to either alone.
The 50% tract and interior surgery deliberately differ from the direct
replication and are reported as transfer tests.

**Matched random controls, 160 pairs.** Each controlled repeat pair has a
random-background counterpart with identical lengths, event direction,
position, gap length, source/target substitution masks and nonzero cyclic base
replacement offsets. Inserted control strings are random with the same length.
The source substitutions leave the control background random. Controls are
reported separately, never pooled to rescue a failed repeat endpoint.

Record exact sequences, motif, effective period, source and target lengths,
all edit geometry, masks/counts, generation seeds and input hashes. The 80
perfect repeat pairs and 80 perfect controls have analytic edit distance equal
to their length difference. Noisy generation edit counts are upper bounds, not
assumed optima.

## Validation and measurement

An independent ordinary unrestricted two-row dynamic program computes exact
edit distance for **all 380 pairs**, before timing. It imports no engine code
and is checked against the Python reference on exhaustive small inputs and
additional random cases. All returned scripts and all finite round intervals
are checked against these independent scores, including warmups. Preserve
budget fallbacks, failures and slow cases; do not silently drop them.

Use unchanged v4 flags: initial width 8, predictor block 64, seed k=11, epsilon
0, order block 8, dense index limit 4,000,000 cells, sparse work limit 4,000,000
expansions, shared bound payload limit 64 MiB, default main DP limit 50,000,000.
Use only the linear predictor. Freeze engine, source, harness, input and
configuration hashes before measurements. Engine binary SHA256 is
`4a492041691bd03ab0577152cfe0aa056a19f8f4be18551e43538d4a8d476a25`.

Rotate/reverse the four setting orders across seven timed repetitions. Shuffle
case order reproducibly per repetition, using the same order for all settings
within that repetition. Preserve one excluded complete warmup per setting.
Internal `total_ms` includes prediction, all bound work, main DP and traceback;
batch wall time is separate. Run no competing correctness jobs while timing.

Report sums of per-case medians, per-case wins/losses, paired ratios and
uncertainty. Bootstrap complete paired clusters, 9,999 resamples, seed 20261009:
20 seed clusters in direct replication, 20 locus clusters in the controlled
suites. All lengths/variants belonging to a cluster move together. Intervals
describe input-cluster variation in this synthetic sample, not uncertainty over
machines or biological populations. Do not count timing repetitions as new
biological observations.

## Secondary diagnostics and decision

Report controlled-repeat transfer and random controls separately, and
predefined strata by length, motif period, gap fraction, direction and regime.
Include alignment cells, rounds, cells after the first optimal upper bound,
preprocessing/certificate times, sparse expansions and cache states, payload
memory, and all fallback counts. Check sparse/dense alignment-work agreement
when budgets allow it. Components and subgroup observations are explanatory;
they do not redefine the primary endpoint.

If the primary gate fails, record failure even if a subgroup is promising. If
direct replication succeeds but transfer fails, describe a restricted synthetic
mechanism rather than a broadly useful repeat aligner. No external production
aligner comparison, affine-gap scoring, real-read or genotype-to-phenotype
benefit is claimed by this experiment.
