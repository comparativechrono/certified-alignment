# Prototype v3: ordered certification

Run date: 9 October 2026.

The ordered block/gap bound is implemented and passes the independent audits.
It reduces expanded alignment work, including repeatable savings on some
compensated-indel cases. **Its dense preprocessing costs more than it saves in
aggregate, and the narrower successes remain slower than ordinary linear
banding.** The default remains the original length bound.

This experiment does not establish a useful performance advantage for the
mathematics-inspired cost predictor. The conclusion concerns this explicit
surrogate and these implementations; it does not test the original theoretical
construction in full or other possible transfers to bioinformatics.

## Implementation

`--bound ordered` takes the maximum of the v2 histogram bounds and an ordered
block relaxation. Consecutive source blocks align to consecutive target chunks.
An exact block match costs zero; other chunks pay at least one and at least their
length imbalance. The relaxation allows every target placement, including
repeated matches, without forcing any predicted anchor.

An independent proof gives a three-candidate recurrence and constant-time queries
at arbitrary suffix positions, including inside a source block. See
`ORDERED_BOUND.md`. Fixed-size source blocks of eight bytes give roughly one
index row per eight source positions. Exact packed word identities avoid
collision assumptions and repeated-hit enumeration.

The index is built lazily when the cheaper bounds are insufficient, once per
pair, and retained across expansion rounds. It has O(nm/k) preprocessing and
storage: quadratic for the fixed k=8 used here. Default limits are 4,000,000
ordered cells and 64MiB shared bound-index payload. Cap exhaustion safely omits
this component and reports the fallback. All preprocessing is timed.

The predictor source is byte-for-byte unchanged. The seed/gap reasoning has
prior art in A*PA, cited in `ORDERED_BOUND.md`; this is an experimental combination
with our band engine, not a claim of a newly discovered alignment theorem.

## Independent correctness checks

All checks passed:

- 28,665 end-to-end output records checked against independent unrestricted
  edit-distance DP: 51,604 finite certificate intervals and 28,661 CIGAR replays.
  There were 28,003 exact outputs, 658 valid explicitly requested approximations,
  and four expected budget errors.
- 16,980 deterministic comparisons against untouched v2 confirmed existing
  modes retain their results. The ordered bound never decreased a matched qgram
  certificate bound or increased the number of alignment rounds.
- 450,184 suffix checks against independent edit distance and a slow relaxed
  oracle that enumerates every target chunk length, rather than using the
  optimized recurrence.
- 483,080 combined-bound checks, 2,250,920 threshold checks, 404,512 Lipschitz
  checks, and 9,568 fallback checks. Cases include all block sizes 1..8, partial
  terminal blocks, arbitrary byte values, empty strings and exact cap boundaries.
- The original 8,813-output alignment test suite also passed.

Reports: `results/ordered_bounds.json`, `results/ordered_integration.json`, and
`results/correctness.json`. Actual allocator failure was not fault-injected;
explicit memory and work-cap fallbacks were tested. The bound's admissibility
comes from the proof; finite tests support its implementation.

## Fixed benchmark

The 148 cases comprise the unchanged 118 inputs from v2, plus 30 new compensated-
indel pairs fixed before inspecting v3 timings. The new cases have source length
2,048, seeds 211/223/227/229/233, event lengths 32/128/384, and either random or
repeated-tract backgrounds. The prior 118 are known test cases at this stage;
only these additional 30 are new to the v3 experiment.

Four unchanged predictors are crossed with length, qgram and ordered bounds.
Settings are fixed at width 8, predictor block 64, seed-k 11, order-block 8,
order-max-cells 4,000,000 and exact epsilon 0. Each setting has one warmup and five
timed repetitions, with rotated/reversed setting order. There was no parameter
tuning after seeing the results.

All 8,880 timed outputs and 1,776 warmup outputs passed certificate, CIGAR and
cross-setting score checks. Forty prior pure-indel constructions also have an
analytically known optimum. The benchmark harness does not run an independent
full-DP oracle on the other pairs; those oracle checks are in the separate audit.

Times sum per-case medians over five repetitions. Total time includes prediction,
all bound preprocessing, certificate evaluation and traceback, but excludes
process startup and JSON output. Small timing differences in this shared runtime
are not statistical significance claims. Main DP cell counts are deterministic.
All comparisons below use the same v3 executable.

## Aggregate result: stronger bounds, excessive index work

Ordered versus qgram on the prior 118 pairs, holding the predictor fixed:

| Predictor | Main DP cell change | DP + pilot + index cell change | Observed total time change |
| --- | ---: | ---: | ---: |
| Linear | -6.86% | +59.61% | +24.87% |
| Cost | -6.22% | +54.83% | +20.65% |
| Anchor-linear | -9.35% | +43.96% | +15.24% |
| Anchor-cost | -8.01% | +44.70% | +12.27% |

The combined cell count is a partial work proxy. Index, pilot and main DP cells
have different execution costs, and this sum omits seed lookup, histogram and
certificate operations. End-to-end time includes those operations.

For anchor-cost, the prior-suite main DP reduction is 125,084,473 to 115,070,979
cells. It costs another 68,192,403 ordered-index cells. Measured ordered
preprocessing totals 92.708ms, while measured DP time falls by 40.626ms. With
other overhead included, total time rises from 497.908 to 559.009ms.

On the 30 new pairs, anchor-cost saves only 1.65% of main DP cells (78,794,495
to 77,490,984), adds 15,797,790 index cells, and total time rises from 278.019
to 291.366ms, or 4.80%. Larger and repetitive compensated indels remain difficult.

No index-cap fallback occurred in the benchmark. The largest ordered payload
was 9,494,572 bytes, excluding allocator bookkeeping and other process memory.
Fallback behavior therefore does not explain the negative aggregate result.

## The narrower positive result needs the right baseline

For all five prior 64-base compensated-indel cases, anchor-cost with ordered
certification roughly halves main DP cells and reduces observed runtime by
about 19% relative to the same predictor with qgram certification. The five
new random-background 32-base cases also roughly halve main DP cells.

However, plain linear banding is already faster in both groups:

| Five-case group | Linear / length (ms) | Anchor-cost / qgram (ms) | Anchor-cost / ordered (ms) |
| --- | ---: | ---: | ---: |
| Prior compensated, 64-base events | 8.225 | 19.537 | 15.825 |
| New random compensated, 32-base events | 3.921 | 11.624 | 11.260 |

For the 64-base group, linear uses 2,471,980 main DP cells. Ordered anchor-cost
uses 2,483,352 main cells, plus 184,255 pilot and 2,632,965 ordered-index cells.
For the 32-base group, the corresponding counts are 1,171,175 main cells for
linear, and 1,179,293 main + 169,388 pilot + 2,632,965 index for ordered anchor-cost.

The new certificate removes an extra expansion induced by the anchored band.
Plain linear already has zero post-optimum expansion in these groups. Ordered
anchor-linear takes 16.438ms and 11.095ms respectively, so the evidence also
does not isolate a robust benefit from cost prediction over identical-anchor
linear interpolation.

## The difficult score-818 example remains unresolved

On `compensating_indels-n4096-s17`, anchor-cost finds the final score 818 at radius
128. The ordered lower bound improves from the qgram value 42 to 268 at that
radius, and from 368 to 477 at radius 512. It still cannot certify before
radius 1,024. Main DP remains 15,309,953 cells, with another 2,101,761 index cells.

This fits a known limitation of the relaxation: on equal-length strings, pairing
every block with an equal-sized target chunk costs at most ceil(n/k), here 512.
The method preserves order but undercharges heavily edited blocks. More
information in the bound alone is not enough to justify its preprocessing cost.

## Decision and possible next experiment

Keep this bound as an experimental option and retain the simple default. Further
polishing of the dense index is not the priority supported by these results.
The original cost surrogate has not demonstrated a useful advantage over linear
banding across the three experiments.

If this direction continues, the next distinct hypothesis should be a sparse or
incrementally constructed order-aware certificate that avoids most of the dense
preprocessing, with ordinary linear banding retained as an essential baseline.
Only after an end-to-end gain against that baseline should we invest in comparison
with optimized external aligners and real sequencing data. Such a sparse version
has not been implemented or tested here.

## Reproduce

Tested on Linux x86_64 with GCC 13.3.0, Python 3.12.14, and
`-O3 -DNDEBUG -std=c++17 -Wall -Wextra`. Python dependencies are standard-library
only. From the package root:

```sh
make
make test
make test-ordered
make benchmark-ordered
./build/certalign --mode anchor-cost --bound ordered < examples/pairs.tsv
```

Reports regenerate in place. Preserve a copy of supplied results before rerunning
if you need the recorded measurements. `reference_v1/` and `reference_v2/` contain
unchanged earlier source for regression checks. Platform-specific binaries are
not packaged.

`results/ordered_benchmark/` includes every input, generation metadata, raw JSON
and CSV outputs, round traces, fallback counts and matched summaries. Earlier
results are retained separately. `MANIFEST.sha256` records the packaged files.

Measured executable SHA256:
`f88a3d9c7c3b5e746fc15d07cc09979724aa73f74e57d69dad45898389d5eabf`.
Combined input SHA256:
`80e6c6e100848eb13a364a9395fec7e82191cbbe69ca8d47444c619b78249579`.
Predictor SHA256, unchanged from v1:
`930a1a18658e3a0b8b651291199afb93912054672daafd7163eb20a112d1156a`.

All experiments use synthetic pairs and global unit-cost scoring. No external
aligner, affine-gap model, real sequencing dataset or genotype-to-phenotype
association method was evaluated.
