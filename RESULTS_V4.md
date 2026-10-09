# Sparse threshold certificate: v4 results

Run date: 9 October 2026. **The speed objective failed.** The sparse evaluator
preserved the dense ordered certificate's stopping rounds and used much less
bound memory, but repeated threshold traversal cost more than it saved. No
sparse setting beat simple linear banding in aggregate on either the previous
148 pairs or the 30 new indel pairs. These are results for this implementation
and these synthetic inputs, not an impossibility result about sparse proofs.

## Controlled change

V4 changes only evaluation of the ordered relaxation: lazy threshold proofs,
cached lower/upper intervals, a cheap feasible upper-bound guard, and an explicit
stack replace dense score-table construction. `--bound sparse` retains all q-gram
bounds. Predictors and the main first-exit certificate are unchanged. The dense
`ordered` mode remains available. See `SPARSE_BOUND.md` for the proof.

The cache stores only established inequalities. Budget exhaustion is unknown;
it cannot certify an unfinished proof or establish a new upper bound. Allocated
tiles have 64 columns. The directory is still partly dense, and worst-case work
is not newly subquadratic. Revisiting a state consumes work each time.

## Fixed experiment

- 178 pairs: all 148 v3 pairs unchanged, plus 30 new compensated-indel pairs.
- New inputs have source length 3,073, gap lengths 47, 193 and 511, random and
  repeated backgrounds, and seeds 401, 409, 419, 421 and 431. They were generated
  and frozen before v4 timing. The old 148 are known data, not a fresh test set.
- Four unchanged predictors crossed with length, qgram, ordered and sparse:
  16 settings, one complete warmup and five timed repetitions. Setting order
  rotates by five positions and reverses on odd repetitions.
- Initial radius 8, predictor block 64, seed k=11, exact epsilon=0, order block
  8, dense limit 4,000,000 cells, sparse limit 4,000,000 expansions, shared bound
  payload cap 64 MiB. No tuning after observing results.
- GCC 13.3, `-O3 -DNDEBUG -std=c++17 -Wall -Wextra`, Linux x86_64, Python 3.12.
- Timings below sum each case's median internal `total_ms`. They include
  prediction, preprocessing, proof traversal, allocation, alignment DP and
  traceback. Serialization/process startup are excluded; batch wall times are
  separately retained. Component medians need not sum exactly to total medians.

All 17,088 benchmark outputs passed validation (14,240 timed, 2,848 warmup).
Every script replayed and all settings agreed on scores and exact certificates.
Forty preserved cases have independent analytic optimum scores. The benchmark
does not run a full-DP oracle for the other long inputs; independent oracle
coverage is described below. Total benchmark wall time was 116.79 seconds.

## Primary comparison: ordinary linear banding

The fastest sparse setting over both main scopes was **linear/sparse**. This is
a descriptive selection among four fixed variants using observed aggregate
timings, not a separately validated selection rule. No per-pair winner oracle
is used.

| Scope | Linear/length (ms) | Linear/sparse (ms) | Runtime change | Main DP cell change |
|---|---:|---:|---:|---:|
| Previous 148 | 424.904 | 614.016 | +44.5% | -9.27% |
| New 30 indel pairs | 210.825 | 265.042 | +25.7% | -2.06% |
| All 178 | 635.729 | 879.059 | +38.3% | -6.83% |

Linear/sparse was slower on 137 of the previous 148 individual case medians and
29 of the new 30. Small individual timing differences should not be overread
in this shared runtime. On all 178, main DP cells fell from 216,216,783 to
201,459,363, but sparse certificate time was 281.591 ms versus 12.249 ms for
linear/length. Linear/qgram took 635.806 ms overall: effectively the same
observed aggregate as length banding, with no established runtime gain.

The cost-informed anchored variant was further behind: anchor-cost/sparse took
1,227.312 ms on the previous 148 and 845.477 ms on the new 30, versus the
424.904 ms and 210.825 ms linear/length baselines. This experiment supplies no
evidence that the cost predictor improves this approach.

One narrower signal is worth preserving: on the **nine known tandem-repeat
cases**, linear/sparse took 7.932 ms versus 18.971 ms for linear/length, 10.287 ms
for linear/qgram, and 13.807 ms for linear/ordered. Main DP cells fell from
6,525,969 for linear/length to 1,539,228. This subgroup was inspected after the
aggregate comparison and comes from the original dataset; it is not fresh
confirmation. None of the new v4 gap/background groups beat linear/length.
The new repeated-background compensated-indel cases are a different generation
model, so they neither validate nor precisely replicate the old tandem-repeat
signal. A new repeat-specific dataset would be needed to pursue that lead.

## Did sparsity improve the dense bound?

For every one of the 712 case/predictor combinations, sparse and dense ordered
evaluation used exactly the same cumulative alignment DP cells. There were no
work-cap, memory-cap or allocation fallbacks in either mode. Thus weaker
certification from budget exhaustion does not explain these timings.

| Predictor | Dense ordered, previous 148 (ms) | Sparse, previous 148 (ms) | Change |
|---|---:|---:|---:|
| Linear | 527.080 | 614.016 | +16.5% |
| Cost | 558.719 | 627.105 | +12.2% |
| Anchor-linear | 798.093 | 1,181.186 | +48.0% |
| Anchor-cost | 789.757 | 1,227.312 | +55.4% |

On the new 30, linear/sparse was 1.3% faster than linear/ordered (265.042 versus
268.598 ms), a small observed difference. It still lost to linear/length by
25.7%. Anchor-cost/sparse was 48.4% slower than its dense counterpart.

The anchor-cost component breakdown shows where the savings went:

| Component | Dense, previous 148 (ms) | Sparse, previous 148 (ms) | Dense, new 30 (ms) | Sparse, new 30 (ms) |
|---|---:|---:|---:|---:|
| Bound preprocessing | 112.563 | 4.038 | 44.597 | 1.443 |
| Certificate | 38.216 | 545.730 | 21.788 | 312.895 |
| Alignment DP and traceback | 571.251 | 606.332 | 482.876 | 505.138 |
| Total engine time | 789.757 | 1,227.312 | 569.681 | 845.477 |

Anchor-cost sparse expanded 31,321,779 proof nodes on the previous 148 while
caching 1,473,823 distinct states: about 21 expansions per cached state in
aggregate. On the new 30 it expanded 17,590,470 nodes and cached 964,591 states.
There were 60,018,207 and 33,558,681 cached decisions, respectively. Threshold
revisits and their control flow are substantial despite the sparse footprint.
Sparse expansions and dense recurrence cells have different costs; adding them
is only a partial work proxy, not a runtime model.

The cheap upper guard rejected 2,975 and 3,040 node-threshold requests in these
two anchor-cost scopes. These are ordinary disproofs, not limit fallbacks.
Maximum per-pair sparse expansions were 3,982,251, below the fixed 4,000,000 cap.

## The genuine benefit is bound memory

These are maximum allocated **extra ordered/sparse payloads per pair** in each
scope, not process RSS and not sums across cases. Sparse bytes include its
directory, exact packed words, full-depth stack and tiles. The main alignment
score matrix and allocator overhead are excluded from both figures.

| Anchor-cost scope | Dense ordered bytes | Sparse bytes | Reduction |
|---|---:|---:|---:|
| Previous 148 | 9,494,572 | 1,524,144 | 83.9% |
| New 30 | 4,795,464 | 1,021,352 | 78.7% |

Including histogram payload, peak bound bytes were 9,539,660 versus 1,567,696 on
the previous 148, and 4,832,878 versus 1,058,766 on the new 30. This is a useful
engineering property, but it has not demonstrated an overall memory advantage
against simple length banding, nor a useful speed/memory tradeoff on real data.

## Independent correctness audit

The new suffix test uses a slow minimum over **all possible target chunk
lengths**, rather than copying the optimized three-candidate recurrence.
It checks thresholds and actual cached interval endpoints after successful,
failed and interrupted proofs. The completed checks were:

| Check | Count |
|---|---:|
| Independent suffix oracle checks | 302,704 |
| Uncapped threshold checks | 1,816,224 |
| Capped threshold checks | 272,676 |
| Cached lower/upper endpoint checks | 116,999 |
| Complete integration outputs | 41,035 |
| Finite certificate intervals | 73,121 |
| CIGAR replays | 41,031 |
| Unchanged-v3 deterministic comparisons | 22,640 |
| Uncapped sparse/dense exact comparisons | 8,530 |

All passed. Cases cover block sizes 1..8, arbitrary bytes, empty strings,
partial source blocks, tile boundaries, shuffled thresholds, repeated queries,
work and memory limits, and checked-size overflow rejection. Nonexact
approximation results (608) retained their requested guarantees. Four expected
main-DP-budget errors remained explicitly uncertified.

The original `make test` was also rerun: 8,813 outputs and 23,580 finite
certificate intervals passed against unrestricted edit distance. Actual
allocator-failure injection was not performed; memory-cap interruption was
tested independently. The implementation and proof received a separate manual
review with no material defect found.

## Research decision and limits

Keep the sparse implementation as a tested, lower-memory reference, with the
ordinary length bound still the default. Do not present v4 as a performance
improvement or a newly enabled bioinformatics algorithm. This controlled test
rejects the current hypothesis that sparse threshold evaluation alone rescues
the ordered bound's runtime.

Further work would need a distinct, falsifiable mechanism: for example, sharing
work across changing thresholds or using a stronger inexpensive indel bound.
It should first demonstrate an advantage against linear/length on fresh data
before adding biological complexity. Continuing to tune the current branch on
these known pairs would not provide that evidence. The known tandem-repeat
subgroup provides a narrower candidate for fresh validation, rather than a
reason to claim general indel-alignment improvement.

All experiments are global unit-cost synthetic alignment. No affine gaps,
real-read benchmark, external production aligner comparison, or downstream
genotype-to-phenotype benefit has been tested. The mathematical release inspired
the original surrogate predictor; these experiments do not establish that its
theorem transfers to this prototype. Provenance remains in README.md.

## Reproduction and frozen identities

```sh
make
make test
make test-sparse
make benchmark-sparse
```

`reference_v3/` preserves the previous engine source for regression checks.
`results/sparse_benchmark/` contains the frozen input configuration, exact pairs,
generation metadata, all raw outputs, component timings, trace diagnostics and
paired summaries. Input/configuration files refuse silent modification on a
rerun. New timings can differ with the machine and runtime environment.

- Engine binary SHA256: `4a492041691bd03ab0577152cfe0aa056a19f8f4be18551e43538d4a8d476a25`
- Unchanged predictor SHA256: `930a1a18658e3a0b8b651291199afb93912054672daafd7163eb20a112d1156a`
- All 178 input pairs SHA256: `96492c28420e4eb6636f46863478676976a227431cf7eb11ceeed073bdb4ceb6`
- Frozen configuration SHA256: `b61fe4077bfbcfa431721872a8580f46c253fb5d8b79c525be6ba3038fe7f4b7`

Binaries are excluded from the source archive; these hashes identify the run.
