# Prototype v1: measured results

Run date: 9 October 2026.

The prototype implements a cost-informed alignment band with an independent
exactness certificate. Correctness checks passed. **This first cost surrogate
did not establish a runtime improvement over matched linear-band baselines.**
The certificate and experiment harness are usable foundations for further work;
the current predictor is not ready to recommend as a faster alignment method.

## Scope and hypothesis

The experiment asks whether a finite, local edit-cost proxy can place a global
unit-cost alignment band more efficiently than linear interpolation. It tests
both endpoint-only proposals and proposals using identical unique-k-mer anchors.
Every benchmark run requests an exact result (`epsilon=0`).

This is an experimental adaptation inspired by the mathematics linked in the
README. It does not implement the paper's implicit minimum-prefix-cost predictor
or inherit its runtime guarantee. The exactness guarantee comes from the separate
first-exit certificate whose proof appears in the README.

## Correctness evidence

The independent Python full-DP oracle checked 8,813 output records:

- 961 exhaustive pairs over binary strings of length 0 through 4, plus 171
  random/adversarial pairs, exercised across modes and parameter settings.
- All 23,580 finite intermediate certificate intervals contained the oracle
  optimum. There were also 1,321 rounds without a finite within-band path.
- All 8,808 successful outputs had valid replayed edit scripts: 7,986 were
  certified exact and 822 were nonexact results meeting their explicitly
  requested approximation guarantee.
- Five deliberate budget-exhaustion cases returned errors, not false
  certification. There were no test failures.

The full record is `results/correctness.json`. Testing supports the implementation;
the mathematical guarantee is the certificate proof, subject to correct code.

## Benchmark design

There are 63 synthetic pairs: seven scenarios, source lengths 256, 1,024 and
4,096, and seeds 17, 29 and 43. Scenarios cover uniform errors, clustered errors,
long insertions, long deletions, tandem repeats, compensating indels, and high
similarity. Each mode was run three times. Initial radius was 8, block length
64, and seed length 11. The dataset and generator are included.

All four banded variants ran on all 63 pairs. Full DP ran on the 33 pairs for
which both sequences had length at most 1,024. There were 855 measured output
records with no validation failures. Every result was certified exact; edit
scripts and scores agreed across modes. The other 30 pairs have certificate
and cross-mode evidence, not an independent full-DP reference in this benchmark.

For each case and metric, take the median of three runs. The table sums those
per-case medians across the same 63 cases. Time is internal total wall time,
including prediction and traceback but excluding process startup and JSON
serialization. DP cells include repeated work when bands expand. Pilot cells
are additional prediction work and must not be omitted when assessing cost.

| Proposal | Main DP cells | Additional pilot cells | Total time (ms) | Total rounds |
| --- | ---: | ---: | ---: | ---: |
| Linear | 71,782,254 | 0 | 230.661 | 209 |
| Cost | 71,806,981 | 5,788,288 | 253.329 | 209 |
| Anchor-linear | 87,098,406 | 0 | 323.541 | 215 |
| Anchor-cost | 84,833,977 | 3,061,779 | 324.741 | 213 |

Interpretation:

- **Cost versus linear:** main DP work increased by about 0.034%, with identical
  total expansion rounds. Adding pilot cells increases the combined cell count
  by about 8.1%. Observed runtime was 9.8% higher. There is no evidence of a
  benefit from this predictor on the endpoint-only comparison.
- **Anchor-cost versus anchor-linear:** main DP cells fell by 2.6%, but adding
  pilot work makes the combined cell count about 0.9% higher. Observed runtime
  was 0.4% higher, effectively tied at this measurement precision. The modest
  reduction in main alignment work did not produce an end-to-end gain.
- Anchors themselves did not improve the aggregate result in this prototype.
  Comparing anchor-cost only with linear would mix two changes and obscure
  the contribution of cost prediction.

Pilot recurrence cells and main allocated DP cells have different overheads,
so their sum is a useful work indicator, not a hardware-independent runtime
model. These short, shared-environment measurements and three repetitions do
not support a statistical significance claim or precise small timing differences.
Deterministic cell counts provide the stronger evidence here.

## What this tells us

The current proxy distributes an interval's *net* length imbalance according
to local edit cost. Local cost can reflect substitutions, repeat ambiguity, or
mispaired blocks as well as actual indel displacement. In an equal-length
interval it predicts no displacement at all. This is a specific limitation of
the implemented surrogate, not a negative result for every prefix-cost predictor.

The certificate deliberately uses only a remaining-length lower bound. A band
can already contain an optimal path yet still need widening to exclude a cheaper
path outside it. Better band prediction alone therefore need not reduce exact
alignment work. Retrospective analysis of the first repetition separates these
events (the traces and cell counts are deterministic across repetitions):

| Proposal | Cases finding the final optimum before certification | Main DP cells spent in later rounds | Share of main DP cells |
| --- | ---: | ---: | ---: |
| Linear | 20 / 63 | 18,284,173 | 25.5% |
| Cost | 20 / 63 | 17,475,790 | 24.3% |
| Anchor-linear | 26 / 63 | 56,712,124 | 65.1% |
| Anchor-cost | 26 / 63 | 58,512,910 | 69.0% |

These later rounds are needed by the current certificate. They cannot safely
be skipped just because the final optimum is known retrospectively. For example,
on `compensating_indels-n4096-s17`, anchor-cost already finds the final score 818
at radius 128, but the lower bound is 0. At radii 256 and 512 the lower bounds are
89 and 368. Certification requires radius 1,024. This illustrates weak bounds,
not an expensive boundary scan: computing the certificate itself accounts for
only a few milliseconds in the aggregate; its stopping rule triggers large DPs.

Cost prediction changes no final round count in the endpoint-only comparison.
With anchors it saves one round on just two cases: `long_deletion-n1024-s43`
and `long_insertion-n4096-s17`. Their main DP cells fall from 492,733 to 246,064
and from 4,066,494 to 2,031,460, respectively; the added pilot work is 16,384
and 49,152 cells. These are substantial isolated savings, each on one seed,
and deserve investigation. They do not establish a general benefit on indels.

The next focused experiment should strengthen the admissible suffix bound while
holding these proposals fixed, then revisit local displacement prediction.
Changing both together would make their effects harder to identify. Preserve
the exact certificate as an independent correctness layer. Any stronger bound
must be proved admissible before being used to certify.

## Reproduce and inspect

Tested with GCC 13.3.0, `-O3 -DNDEBUG -std=c++17 -Wall -Wextra`, Python 3.12.14,
Linux x86_64. No third-party Python packages were used. From the package root:

```sh
make
make test
make benchmark
```

These commands regenerate results in place. The archive retains the original
results for comparison; make a copy before rerunning if needed. Source is
included; compilation is required. Timings depend on the host and compiler.

Useful files:

- `README.md`: interface, algorithm, certificate proof, limitations and provenance.
- `include/predictor.hpp`: experimental predictor, separable from the engine.
- `src/aligner.cpp`: exact band DP, boundary certificate and traceback.
- `tests/test_correctness.py`: independent oracle and adversarial checks.
- `benchmarks/run_bench.py`: deterministic input generation and comparisons.
- `results/benchmark/pairs.tsv` and `cases.json`: all input pairs and metadata.
- `results/benchmark/raw_results.json`: outputs including per-round traces.
- `results/benchmark/raw_results.csv`: tabular measurements.
- `results/benchmark/summary.json`: per-case, aggregate and matched comparisons.
- `MANIFEST.sha256`: checksums of the packaged source and results.

The benchmark records SHA-256 fingerprints for the tested executable and input
pairs. The executable fingerprint was
`58e51f883fd6ccc800b92064100fc2ac107e522f3fe8807c804ca328a2d45eeb`;
the input-pairs fingerprint was
`2ef3d5b6c677b87a81500d778a5e0352c2e5c9d1ac3d52c5cd3f43b620ede585`.

No claims are made about superiority to Edlib, A*PA2, QuickEd, WFA2 or BiWFA:
they were not benchmarked. No real sequencing dataset, affine-gap scoring,
local alignment, or genotype-to-phenotype association method was tested.
