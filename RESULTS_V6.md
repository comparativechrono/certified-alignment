# v6: external comparison of the unchanged sparse-certificate prototype

## Decision

**The practical performance gate failed decisively.** On the primary isolated-indel
workload, the sparse prototype took 29.56 times as long as Edlib and 116.76–125.78
times as long as WFA2. All returned scores, paths and finite certificate intervals
passed validation. The result supports retaining the certificate as research
infrastructure, while stopping performance development of this implementation
as a competitor to these exact unit-cost aligners.

This evaluates the fixed **linear band predictor with sparse certification**.
It does not evaluate every cost-based predictor, affine-gap scoring or a new
alignment engine. The engine and predictors are byte-for-byte identical to v5.

## Primary result

Times below cover all 192 primary pairs. They sum each pair's median complete
native alignment time across seven repetitions and add median batch setup and
teardown. A smaller number is better. Every method returns a full alignment path
under identical global unit-cost scoring; none receives the oracle score.

| Method | Genomic-background total (ms) | Matched random-control total (ms) |
| --- | ---: | ---: |
| WFA2 high memory | 55.153 | 56.365 |
| WFA2 ultralow memory | 59.412 | 60.398 |
| Edlib | 234.653 | 236.896 |
| Prototype sparse | 6,937.078 | 6,839.281 |
| Prototype q-gram | 7,423.546 | 7,414.518 |
| Prototype length | 7,526.983 | 7,423.173 |

The primary endpoint required at least a 10% reduction in runtime against
**each** external configuration, with each paired locus-bootstrap 95% upper
ratio below one and no correctness or budget failures. No external comparison
met the speed criterion:

| External baseline | Sparse / baseline runtime | Paired locus-bootstrap 95% interval |
| --- | ---: | ---: |
| Edlib | 29.56 | 27.92–31.25 |
| WFA2 high memory | 125.78 | 118.85–132.93 |
| WFA2 ultralow memory | 116.76 | 111.09–122.47 |

The actual observed batch totals tell the same story. Their medians were
7,087.075 ms for sparse, 239.560 ms for Edlib, 61.450 ms for WFA2 high and
63.887 ms for WFA2 ultralow. The main conclusion does not depend on summing
per-pair medians instead of taking a median of whole-batch totals.

## What the internal improvement means

Sparse reduced runtime by 6.55% against the prototype's q-gram bound and 7.84%
against its length bound on the primary suite. These are diagnostic internal
comparisons, not the external success criterion. The matched random controls
showed a 7.76% reduction against q-gram, while remaining 28.87 times as slow as
Edlib and 113.24–121.34 times as slow as WFA2.

Sparse cut primary main-DP work by 42.46% relative to q-gram
(1,623,182,212 to 933,945,038 cells). However, its summed median certification
time was 2,935.027 ms, compared with 128.581 ms for q-gram. The stronger bound
saves DP work but is expensive to evaluate; reduced cell count is not equivalent
to a competitive complete alignment algorithm.

Fourteen primary cases and thirteen controls reached the sparse proof-work cap
and used the existing admissible fallback. They still returned exact alignments.
No case was removed, rescued with larger limits, or credited as a speedup for
ending in an error.

Peak whole-process RSS on the primary suite was 118,632 KiB for sparse, 8,448 KiB
for Edlib, 10,960 KiB for WFA2 high and 8,448 KiB for WFA2 ultralow. These measurements
include preloaded inputs, retained paths/diagnostics and library workspaces;
they are not per-alignment auxiliary-memory measurements. This experiment shows
no process-memory advantage for the sparse prototype either.

## Workload, measurement and checks

The primary suite uses eight nonoverlapping E. coli K-12 MG1655 reference windows,
centered lengths 1,024/4,096/8,192, an interior insertion or deletion of 2% or 8%
of source length, and two fixed error regimes. A second 192-pair suite replaces
the reference background with random DNA while preserving edit masks and nested
geometry. These are **semisynthetic sequences**, not sequenced reads or a
platform-calibrated error model. The input auditor independently reconstructed
all 384 pairs.

The protocol, source, dependencies, inputs, oracle scores and runtime freeze were
[committed before timing](https://github.com/comparativechrono/certified-alignment/commit/efe18442a98b1296442b773b87e01c884db1f047).
Edlib and both WFA2 configurations use their official native APIs at pinned
commits. Native timings include prediction, preprocessing, certification,
traceback, CIGAR conversion and required cleanup; WFA2 workspace reuse is allowed
within each 192-pair batch and its setup/teardown are counted. Input parsing,
output serialization and independent validation are outside all native timings.

One excluded warmup and seven measured repetitions produced 96 batches and
18,432 alignment outputs. All passed independent score and CIGAR checks against
384 unrestricted-DP oracle scores. Validation also checked 18,432 final lower
bounds, 48,960 round lower bounds and 48,960 round upper bounds. There were zero
validation failures. Separate adapter tests had already passed 19,458 score/path
checks and 3,243 comparisons with the original CLI.

The bootstrap resamples the eight whole loci, retaining all nested lengths and
edit variants together, for 9,999 resamples with a fixed seed. Its intervals
describe locus variation conditional on this machine and workload. They do not
measure uncertainty across machines, genomes, sequencing platforms or scoring
models. Eight reference windows are insufficient to establish general biological
performance, although the observed external runtime gap is large.

## Reproduction and interpretation

See [reproduction commands](EXTERNAL_REPRODUCTION.md), the
[frozen protocol](EXTERNAL_BENCHMARK_PLAN.md), and
[native API settings](benchmarks/EXTERNAL_DEPENDENCIES.md).
[Machine-readable summaries](results/external_benchmark/summary.json),
[validation](results/external_benchmark/validation.json),
[independent results audit](results/external_benchmark/results_audit.json),
[per-pair medians](results/external_benchmark/per_case.json),
[schedule](results/external_benchmark/schedule.json), and
[raw records](results/external_benchmark/raw.jsonl.gz) are retained.
The raw records are losslessly compressed; the restoration script checks both
compressed and original hashes. Original raw SHA256:
`1e0fd6d7676759107cb64ad4461d3784dc4fe8658e69c24cc58a83d498a17824`.

The earlier indel-specific gain justified this external check. This check does
**not** support a practical improvement over established exact aligners, a new
algorithmic capability, or a benefit for genotype–phenotype association studies.
Any future proposed benefit needs a distinct, stated use case and a new frozen
comparison against the appropriate established methods. Further optimization on
these inputs would be development, not independent validation.
