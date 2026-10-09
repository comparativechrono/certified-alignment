# Prototype v2: stronger certificates and an indel-focused test

Run date: 9 October 2026.

The stronger suffix bounds are implemented and validated. They modestly reduce
work on the original benchmark, but **do not reduce alignment cells on any of
the 55 fresh indel-focused cases**. The earlier isolated cost-predictor wins on
indels do not generalize to this extension. This is useful evidence about which
approach to pursue, not evidence of a new competitive alignment algorithm.

## What changed

Added `--bound composition` and `--bound qgram`, with the latter taking the
maximum of remaining-length, symbol-count, 2-gram and 4-gram lower bounds.
Five separate suffix cursor streams maintain word-count differences efficiently
along band exits. Preprocessing is lazy; its time and retained index payload are
reported. A 64MiB index-payload cap can omit a stage safely, weakening the bound.
The first-exit certificate remains the correctness mechanism. See `BOUNDS.md`.

All four predictors, their parameters, the scoring model and the original input
pairs were held fixed. The predictor source is byte-for-byte identical to v1.
`--bound length` remains the default: the experiment does not support promoting
the stronger bounds or the cost predictor as a general performance improvement.

## Independent validation

- 19,107 end-to-end outputs were checked against independent unrestricted Python
  edit-distance DP. All 35,699 finite intermediate certificate intervals
  contained the optimum, and all 19,092 successful edit scripts replayed.
- Successful outputs included 18,462 exact results and 630 results meeting their
  explicitly requested approximation tolerance. Fifteen expected budget errors
  were reported without false certification.
- 5,660 deterministic results under the length bound matched the untouched v1
  executable. Another 5,660 checks confirmed default and explicit length modes
  agree. Timing fields are excluded from these comparisons.
- Direct suffix tests checked 211,098 queries plus 1,055,490 threshold checks
  across 1,044 sequence pairs. They include all suffixes of tiny exhaustive
  binary pairs, backward cursor movement, resets, arbitrary byte alphabets,
  guaranteed hash collisions and partial memory-budget fallbacks.
- The original 8,813-output test suite also passed on the updated engine.

There were no failures. The reports are `results/bounds_integration.json`,
`results/lower_bounds.json`, and `results/correctness.json`. Tests support the
implementation; the proof in `BOUNDS.md` establishes why the bound is admissible.

## Fixed benchmark

There are 118 pairs:

- The exact original 63 pairs, covering seven scenarios at source lengths 256,
  1,024 and 4,096 and seeds 17, 29 and 43.
- 45 fresh pairs: source length 2,048, seeds 101, 103, 107, 109 and 113, and
  insertion/deletion/compensating-indel cases with event lengths 16, 64 and 256.
- 10 repeat controls: those five seeds, with either a 64-base insertion or a
  64-base deletion in a repeated tract.

The extension, parameters and analysis groups were fixed before inspecting the
new benchmark results. It is a synthetic extension, not real sequencing data.
Every case runs under four predictors and three bounds, with one warmup and
five timed repetitions. Setting order is rotated and reversed across repetitions.
No parameters were tuned after seeing the outcomes.

All 7,080 timed and 1,416 warmup outputs passed certificate, CIGAR and score-
agreement checks. Forty fresh cases have an independently known exact score:
the constructed insertion/deletion path meets the remaining-length lower bound.
The other benchmark cases use certificates and cross-setting agreement; this
harness does not run an independent full-DP oracle on them. Independent oracle
testing is the separate validation described above.

Times below sum each case's median internal wall time over five repetitions.
They include predictor, bound preprocessing, certificate and traceback overhead,
but exclude process startup and JSON output. Small percentage differences on
this shared runtime are not statistical significance claims. Deterministic DP
cell counts are the clearer evidence. All comparisons use the same v2 executable.

## Effect of stronger certification

Across all 118 pairs, qgram versus length, with the predictor held fixed:

| Predictor | Change in main DP cells | Change in observed total time |
| --- | ---: | ---: |
| Linear | -3.50% | +2.53% |
| Cost | -3.44% | -2.70% |
| Anchor-linear | -1.19% | +2.73% |
| Anchor-cost | -2.82% | -1.70% |

This aggregate masks the crucial original-versus-fresh distinction. For
anchor-cost specifically:

| Suite | Length-bound DP cells | Qgram-bound DP cells | Length time (ms) | Qgram time (ms) |
| --- | ---: | ---: | ---: | ---: |
| Original 63 | 84,833,977 | 81,207,806 | 310.579 | 301.389 |
| Fresh indel-focused 55 | 43,876,667 | 43,876,667 | 163.605 | 164.724 |

The original-suite reduction is 4.27% in cells and 2.96% in observed time. Its
saved work comes from tandem-repeat and high-similarity cases, not compensating
indels. On the fresh suite, all three bounds give exactly the same main DP cell
count for each of the four predictors. A bound can be tighter without being
tight enough to cross the stopping threshold at one of the doubled radii.

The new bounds are not free. For anchor-cost across all 118 cases, summed
per-case median certificate time increases from 12.158 to 19.676ms, plus
2.529ms of qgram preprocessing. Those costs are already included in total time.
The full summaries retain the composition-only variant as an ablation.

## Does the cost predictor itself help indels?

Compare anchor-cost against anchor-linear using **the same qgram bound**:

| Suite | Anchor-linear DP cells | Anchor-cost DP cells | Change | Observed time change |
| --- | ---: | ---: | ---: | ---: |
| Original 63 | 85,551,714 | 81,207,806 | -5.08% | -3.22% |
| Fresh indel-focused 55 | 42,375,908 | 43,876,667 | +3.54% | -0.10% |

Cost prediction also incurs extra pilot DP: 3,061,779 cells on the original suite
and 1,998,839 cells on the fresh suite. They are not included in the main-DP
columns above; their execution time is included in total time. The near-equal
fresh-suite timing is not a demonstrated benefit despite extra deterministic
work. In the fresh repeat controls, cost prediction uses roughly 127% more main
DP cells than interpolation through identical anchors. The original isolated
long-indel improvements therefore do not establish a robust indel advantage.

## What remains difficult

For `compensating_indels-n4096-s17`, anchor-cost still first finds the final score
818 at radius 128 and still requires radius 1,024 to certify it under every bound:

| Radius | Restricted upper score | Length lower bound | Composition lower bound | Qgram lower bound |
| --- | ---: | ---: | ---: | ---: |
| 128 | 818 | 0 | 17 | 42 |
| 256 | 818 | 89 | 102 | 113 |
| 512 | 818 | 368 | 368 | 368 |

The extra information is valid but far too weak. On the fresh anchor-cost cases,
39.50% of main DP cells are still processed in later rounds after the final
optimal score was first found. This is retrospective identification, not a
license to skip certification. Pure insertion/deletion cases already certify
when their optimal score is found; the remaining delay is in compensating cases.

The next hypothesis should involve **sequence order and displacement**, such as
an admissible bound based on ordered seed matches and the gaps between them.
Word histograms discard much of that information. It is also worth testing
whether selectively skipping the cost predictor in repeat-rich intervals avoids
its observed extra work, but that would be a new predictor change requiring a
new matched test. Neither proposed follow-up has been implemented or validated
here. Any seed-based certificate must allow absent, ambiguous and incorrect
anchors without assuming they belong to an optimum.

## Reproduction and contents

Tested on Linux x86_64 using GCC 13.3.0 and Python 3.12.14, with
`-O3 -DNDEBUG -std=c++17 -Wall -Wextra`. Python dependencies are standard-library
only. From the extracted package root:

```sh
make
make test
make test-bounds
make benchmark-bounds
./build/certalign --mode anchor-cost --bound qgram < examples/pairs.tsv
```

Commands regenerate reports in place; retain a copy of the supplied results if
you want to compare a rerun. `reference_v1/` contains the untouched original
engine source for deterministic regression checks. The archive includes source
and data, not platform-specific binaries.

`results/bound_benchmark/` contains the fixed input pairs, metadata, raw JSON
outputs (including warmups and every round trace), CSV measurements, and grouped
summaries. `results/benchmark/` retains the original v1 benchmark. `RESULTS.md`
is the initial experiment; this document is the subsequent controlled test.
`MANIFEST.sha256` records checksums for all packaged files.

Fingerprints recorded with the measured experiment:

- Predictor: `930a1a18658e3a0b8b651291199afb93912054672daafd7163eb20a112d1156a`.
- Original pairs: `2ef3d5b6c677b87a81500d778a5e0352c2e5c9d1ac3d52c5cd3f43b620ede585`.
- Combined pairs: `f22ff60cf885ebf420da187d86bc33baac3c80dffd0c9e822bd7c6f71b1382b5`.
- Measured executable: `b851d9182316b4ca43325c593c862bc1e404c069fccf3d73cf86ad62ae6ba087`.

No external aligner was benchmarked. These results establish neither superiority
to mature alignment libraries nor performance under affine gaps, local alignment,
real sequencing errors, or pangenome graphs. Genotype-to-phenotype association
methods were not part of this sequence-alignment experiment.
