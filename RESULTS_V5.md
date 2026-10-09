# Fresh tandem-repeat validation: v5

Run date: 9 October 2026. **The fixed primary performance gate was not met.**
The unchanged sparse certificate was substantially faster than length-only
linear banding on fresh replications of the original repeat generator, but its
incremental advantage over q-gram certification was small and uncertain. The
broader controlled-repeat set was slower than q-gram certification. This does
not confirm a useful repeat-specific speed advantage.

V5 changes the evaluation, not the aligner. The v4 binary, all band predictors,
score model and algorithm parameters were preserved. This experiment does not
validate the cost predictor: the previously selected candidate is linear/sparse.

## Fixed decision rule and fresh data

`REPEAT_VALIDATION_PLAN.md` was written before any new engine result or timing.
The candidate had to be at least 10% faster than **both** linear/length and
linear/qgram on the 60-pair direct-replication set, with both paired cluster
bootstrap 95% upper ratio limits below 1, and all correctness checks passing.
Linear/ordered was a secondary comparator. No setting was chosen again on the
new data; no algorithm parameters were tuned.

| Suite | Pairs | Purpose |
|---|---:|---|
| Fresh replication | 60 | Original tandem-repeat generator unchanged; lengths 256, 1,024 and 4,096; new seeds 501..520. |
| Controlled repeats | 160 | Primitive periods 2..6, lengths 1,024/4,096, interior expansions/contractions, approximately 2%/8% copy changes, perfect/imperfect regimes. |
| Matched random controls | 160 | Same lengths, edit positions, gaps and substitution masks as controlled repeats, with random backgrounds. |

All 380 pairs are disjoint from the previous 178 inputs. The controlled suites
have 20 loci and eight variants per locus. Their repeat tracts occupy roughly
half the source, unlike roughly three quarters in direct replication. Perfect
cases have only the constructed indel. Imperfect cases jointly introduce 1%
source-tract interruptions and 0.4% target substitutions; their separate effects
are not identifiable here. Exact masks and geometry are retained in metadata.

An independent audit found a generator-definition mismatch before timing:
the initial filter rejected motifs with any shorter finite-string period,
whereas the written design required motifs that are not exact powers of shorter
words. The filter was corrected to the written design, inputs were refrozen and
independent scores recomputed. Superseded inputs and scores are retained under
`pre_audit_revision/`; no candidate engine outcome informed this correction.

## Primary result

Times sum per-case medians over seven timed repetitions; one complete warmup
per setting is excluded. All figures include the engine's prediction, bound
preprocessing, proof traversal, alignment DP and traceback overhead.

| Fresh replication, 60 pairs | Total time (ms) | Main DP cells |
|---|---:|---:|
| Linear/length | 123.366 | 41,170,531 |
| Linear/qgram | 77.429 | 23,986,724 |
| Linear/ordered | 109.335 | 15,543,989 |
| Linear/sparse | 72.444 | 15,543,989 |

| Sparse comparison | Runtime ratio | Change | Paired cluster 95% interval |
|---|---:|---:|---:|
| Versus length | 0.587 | 41.3% faster | [0.523, 0.650] |
| Versus qgram | 0.936 | 6.4% faster | [0.800, 1.088] |
| Versus dense ordered | 0.663 | 33.7% faster | [0.611, 0.708] |

The q-gram comparison fails both requirements: the observed gain is below 10%,
and its interval includes a slowdown. The result is compatible with a useful
incremental benefit or a modest disadvantage; this sample has not established
the former. Beating length-only banding does not establish that sparse proofs
add enough value beyond the cheaper q-gram bound.

The bootstrap resamples all three lengths within each of 20 seed clusters,
using 9,999 resamples and seed 20261009. It measures synthetic-input variation
conditional on this runtime, not machine-to-machine or biological-population
uncertainty. The two co-primary comparisons were fixed in advance. Subgroup
intervals are descriptive and unadjusted.

| Direct-replication source length | Pairs | Qgram (ms) | Sparse (ms) | Sparse change |
|---|---:|---:|---:|---:|
| 256 | 20 | 0.435 | 0.398 | -8.5% |
| 1,024 | 20 | 3.419 | 2.711 | -20.7% |
| 4,096 | 20 | 73.576 | 69.334 | -5.8% |

The long cases dominate the aggregate. Sparse was faster on nine of the twenty
4,096-base cases and slower on eleven. At length 256, sparse expanded no proof
nodes and used the same DP work as qgram; that very small absolute timing
difference is not evidence of a sparse-proof mechanism. The 1,024-base subgroup
is favorable but does not replace the failed primary endpoint.

## Transfer and controls

These suites were fixed in advance and measured independently of the primary
decision. They are not pooled into a replacement endpoint.

| Suite | Length (ms) | Qgram (ms) | Dense ordered (ms) | Sparse (ms) | Sparse versus qgram |
|---|---:|---:|---:|---:|---:|
| Controlled repeats, 160 | 795.623 | 517.682 | 720.751 | 608.309 | 17.5% slower |
| Matched random controls, 160 | 962.903 | 951.787 | 895.701 | 802.877 | 15.6% faster |

The paired cluster ratio intervals versus qgram are [1.146,1.203] for controlled
repeats and [0.835,0.858] for random controls. Each resamples the same 20 complete
locus clusters, retaining the eight correlated variants. The controlled-repeat
result argues against generalizing the initial repeat-specific gain across
these changes in tract geometry and noise.

The predefined regime breakdown shows a different mechanism worth preserving:

| Background and regime | Pairs | Qgram (ms) | Sparse (ms) | Runtime change | Main DP cell change |
|---|---:|---:|---:|---:|---:|
| Repeat, perfect | 80 | 170.609 | 187.940 | +10.2% | 0% |
| Repeat, imperfect | 80 | 347.073 | 420.370 | +21.1% | -5.2% |
| Random, perfect | 80 | 326.474 | 372.286 | +14.0% | 0% |
| Random, imperfect | 80 | 625.313 | 430.591 | -31.1% | -48.1% |

For the noisy random controls, sparse reduced main DP cells from 207,536,324 to
107,757,368. It eliminated the 99,778,956 DP cells that qgram spent in rounds
after its first optimal upper bound. Sparse and dense ordered both had zero
post-optimal expansion cells on these controls. This directly supports a
certificate-strength explanation for the saved alignment work. It is not
evidence that the cost-based predictor improved.

The noisy random-control ratio interval is [0.678,0.707], but this subgroup is
secondary to the failed repeat hypothesis. Perfect controls instead slow down,
and only 47 of the 80 noisy random case medians were faster despite the weighted
aggregate gain. A future selector cannot assume access to the generator's
perfect/imperfect labels. This is a lead for noisy isolated-indel workloads,
not yet a validated routing rule or production-aligner improvement.

## Overhead, certificates and memory

All 380 case-level sparse/dense comparisons had identical main DP cells, and
neither bound hit a memory or work fallback. Budget weakening therefore does
not explain the outcome.

| Suite | Qgram certificate (ms) | Sparse certificate (ms) | Sparse expansions | Cached states |
|---|---:|---:|---:|---:|
| Fresh replication | 4.822 | 23.520 | 1,092,866 | 444,005 |
| Controlled repeats | 25.126 | 128.899 | 6,329,056 | 1,410,651 |
| Random controls | 26.539 | 171.407 | 9,050,964 | 992,747 |

Sparse proof overhead is substantial. Whether it helps depends on how many
extra band expansions its stronger bound prevents. The table's proof expansions
and DP cells are different operations; their counts are not interchangeable
runtime units. Component medians need not sum exactly to the median total.

Sparse retains its memory advantage over the dense bound. Peak extra payload
across each suite was 867,992 versus 8,614,952 bytes on direct replication,
1,214,616 versus 9,113,472 on controlled repeats, and 1,145,976 versus 9,113,472
on random controls. These are allocated bound payloads, including the sparse
directory/words/stack/tiles, excluding histogram payload, main DP storage and
allocator metadata; they are not process RSS. Maximum per-case sparse work was
618,773 expansions, below the unchanged 4,000,000 cap.

## Correctness and reproducibility

All **12,160 engine outputs** passed: 380 inputs x four bounds x seven timed
repetitions plus one warmup. This includes 10,640 timed outputs. Every distance,
CIGAR and finite round-bound interval was checked against an independent exact
reference. No case, error or slow output was dropped. Benchmark wall time was
58.820 seconds, excluding oracle computation and summary aggregation.

The independent C++ reference imports no engine code and computes unrestricted
two-row Levenshtein distance. It evaluated 3,217,291,264 recurrence cells across
all 380 inputs. Its helper passed 1,168 checks against the separate Python
oracle: 961 exhaustive small binary pairs, 100 random byte pairs, 100 random DNA
pairs and seven edge cases, plus three malformed-input rejection checks. All
160 analytic pure-indel distances matched. The final independent reference pass
took 3.524 seconds and was completed before any algorithm timing.

Setting order rotates and reverses across repetitions. Case order is shuffled
reproducibly within repetitions, identically for all four settings. No competing
validation jobs ran during timing. The original engine and all frozen source
hashes remained unchanged through the run. Environment: GCC 13.3,
`-O3 -DNDEBUG -std=c++17 -Wall -Wextra`, Python 3.12, Linux x86_64.

```sh
make
make benchmark-repeat
```

Reproduction writes to a fresh `results/repeat_reproduction` directory, leaving
the reported data intact. For another run use, for example:

```sh
make benchmark-repeat REPEAT_OUTPUT=results/repeat_reproduction_2
```

The target builds the independent reference, regenerates the fixed cases,
validates all scores, and then runs timings. Portable rebuilds freeze the actual
binary hash in their separate output directory. Existing started experiments
cannot silently be overwritten. The reported run used strict matching to the
unchanged v4 binary; source hashes, seeds, raw records, intervals and diagnostics
are under `results/repeat_benchmark/`. Binaries are not included in the archive.

- Engine SHA256: `4a492041691bd03ab0577152cfe0aa056a19f8f4be18551e43538d4a8d476a25`
- Protocol SHA256: `ebede74f93c727f5a3f7f46eb2a54c294ec3cdea8e4affeedf68015e393fc51b`
- Final pairs SHA256: `193b67af8910b558c78db862da5e22fcff0bce800c481df1a5ece73148fda0f3`
- Final configuration SHA256: `4d3297d40e5555e0373fc286d182dc8f43909b2fb9db86b085cfde5a16933b75`

## Research decision

Do not claim a confirmed repeat-specific algorithmic improvement from this
branch. Preserve the certificate implementation and the noisy-indel control
result, with the failed primary endpoint equally visible. The next useful
question, if pursued, is whether the noisy isolated-indel advantage persists on
independent workloads and against optimized external aligners with the same
scoring model; another selected repeat subgroup is insufficient evidence.

All cases are synthetic global unit-cost alignments. There is no affine-gap,
real-read, external production-library, or genotype-to-phenotype validation.
No newly transferred theorem from the maths release or new asymptotic bound is
established by this experiment.
