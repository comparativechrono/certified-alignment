# Certified cost-based alignment bands

Canonical repository: https://github.com/comparativechrono/certified-alignment.
Future source changes, protocols and results are maintained here. Large
historical raw outputs are stored losslessly as `.gz`; restore them with
`python3 scripts/restore_recorded_results.py` when needed. Reports and inputs
are directly available in the repository.

A reproducible research prototype for **global, unit-cost pairwise alignment**.
Insertions, deletions and substitutions each cost one. The default output is
an exact optimal score and a valid edit script. An optional explicit epsilon
permits a deterministically certified approximation.

The experimental component is a cost-informed proposal for the search band.
Correctness comes from a separate shortest-path certificate, not from the
predictor, anchors, a probabilistic assumption, or a claim that the computed
path avoids a band edge. This package is an initial experiment, not a new
asymptotic alignment theorem or a competitive replacement for mature aligners.

**Current result: v6 external validation failed the performance gate.** The
unchanged linear/sparse prototype was tested on a fresh semisynthetic isolated-indel
comparison. All 18,432 returned alignments across the six methods passed validation,
but the sparse prototype took 29.56 times as long as Edlib
and 116.76–125.78 times as long as WFA2 on the primary workload. Its 6.55% gain
over the internal q-gram baseline does not establish a practical advantage.
See [the v6 report](RESULTS_V6.md), [the protocol](EXTERNAL_BENCHMARK_PLAN.md),
and [reproduction instructions](EXTERNAL_REPRODUCTION.md). The protocol and
inputs were committed before timing; published runs are preserved alongside
source, proofs and earlier results.

**New research: reference-panel adequacy.** The [baseline study](research/PANEL_BASELINE_RESULTS.md)
is complete. On 806 scorable held-out historical leave-one-out calls, the
calibration-selected native filter retained 63 sequence discrepancies at 95%
retention; an ordinary combined confidence model retained 65. This establishes
residual discrepancies, not an advantage for a new algorithm or proof of
missing structural alleles. A small read-level audit is the next feasibility
step. [Frozen protocol](research/PANEL_BASELINE_PROTOCOL.md), source snapshots,
analysis code and independent checks are included.

**Version 4** evaluates that ordered relaxation on demand, proving only the
threshold needed at each band exit and retaining sparse cached intervals.
Work and memory caps preserve admissibility. See `SPARSE_BOUND.md` for the
proof and `RESULTS_V4.md` for the fixed experiment against simple linear
banding. Sparse proof traversal, allocation and preprocessing are included
in total timings.
The fixed v4 experiment found lower bound memory but slower aggregate runtime
than simple linear banding. A faster known tandem-repeat subgroup remains an
unconfirmed lead; no general performance improvement is established.

**Version 5 is a validation release of the unchanged v4 engine.** It tests that
repeat-specific lead on 380 fresh pairs with independent exact distances, a
fixed candidate and a decision rule written before timing. See
`REPEAT_VALIDATION_PLAN.md` for the frozen design and `RESULTS_V5.md` for the
outcome. Source code and scoring parameters are unchanged.
The repeat-specific performance gate failed: sparse was 6.4% faster than qgram
on direct replication with uncertainty spanning a slowdown, and 17.5% slower
on the broader repeat set. A separate noisy random-indel control showed a gain;
it is a secondary lead, not confirmation of the repeat hypothesis.

**Version 3** added an ordered block/gap relaxation with bounded preprocessing,
without changing any predictor. See `ORDERED_BOUND.md` for the proof and
`RESULTS_V3.md` for the new experiment. Its preprocessing cost is included in
all reported total timings.

**Version 2** added independently admissible suffix composition and q-gram bounds
while keeping every band predictor byte-for-byte unchanged. See `BOUNDS.md`
for the inequalities and incremental evaluation, and `RESULTS_V2.md` for the
controlled experiment. `RESULTS.md` records the initial v1 experiment.

## Build and run

Requirements: GCC or Clang with C++17 and unsigned 128-bit integer support;
Python 3.9+ for the tests and benchmark. No Python packages are required.

```sh
make
./build/certalign --mode anchor-cost < examples/pairs.tsv
./build/certalign --mode anchor-cost --bound ordered < examples/pairs.tsv
./build/certalign --mode anchor-cost --bound sparse < examples/pairs.tsv
make test
make test-sparse
make benchmark-sparse
make benchmark-repeat
```

The earlier experiment targets `test-bounds`, `test-ordered`, `benchmark`,
`benchmark-bounds`, and `benchmark-ordered` remain available.

`benchmark-repeat` compiles an independent full-DP reference, checks every fresh
pair, then measures four fixed bounds. It writes to `results/repeat_reproduction`
to preserve the supplied experiment. Each additional run needs a fresh directory,
for example `make benchmark-repeat REPEAT_OUTPUT=results/repeat_reproduction_2`.
The portable reproduction flag permits compiler-dependent binary hashes while
freezing the actual binary and source hashes before that run. The published
v5 experiment used the exact v4 binary.

Input is TSV with exactly three fields: identifier, first sequence, second
sequence. Empty sequences are supported. Matching is case sensitive and treats
all characters literally: N is not a wildcard. There is no reverse-complement
search, local/semiglobal mode, affine gap model, FASTA parser, or sequence graph.

Output is JSONL, one object per input record. `error` means the operation failed
or hit a budget; never interpret that record as an exact result. Check
`certified_exact`, not merely that a score was returned. An error record can
retain the previous round's diagnostic bounds without an edit script.

```sh
./build/certalign --mode cost --width 8 --block 64 --epsilon 0.05 < examples/pairs.tsv
```

`epsilon=0` (the default) requires an exact optimum. The optional approximation
test uses exact rational arithmetic on the requested decimal, not a floating
point tolerance. `distance` is the returned upper bound U, not necessarily the
true optimum when a positive epsilon is requested. CIGAR uses `=` for matches,
`X` for substitutions, `I` to consume a character of the second sequence, and
`D` to consume one of the first. Multiple optimal scripts can exist.

## What is implemented

| Mode | Band proposal |
| --- | --- |
| `full` | Entire DP matrix; an in-process reference, not an optimized library. |
| `linear` | Straight line joining the sequence endpoints. |
| `cost` | Cumulative block-edit-cost proposal between the endpoints. |
| `anchor-linear` | Linear interpolation between sparse unique-k-mer anchors. |
| `anchor-cost` | Cost-informed interpolation between exactly the same anchors. |

Anchors are unique exact k-mers in both strings, put in a monotone chain by
longest increasing subsequence and thinned to spacing `max(32,4*block)`.
They shape the proposal only. They are never forced matches or restrictions
on the global optimum. The certificate can require expanding around wrong
anchors until paths avoiding them are included.

Within each proposed interval, partition the first string into blocks of at
most `block` characters and pair them with proportional chunks of the second.
Compute each block pair's exact edit distance c_j using a small DP. Set positive
proxy mass w_j=c_j+1, with the fixed +1 providing smoothing. For an insertion
imbalance Delta>0, allocate displacement Delta*w_j/sum(w). For a deletion
imbalance, allocate min(block_length_j,lambda*w_j) deletions, with lambda found
by bisection so the required net deletion is allocated. This capacity rule
keeps the center nondecreasing. Interpolate linearly within each block.

This is a **finite explicit surrogate**, inspired by the prefix-cost predictor
in family 121. It is not the paper's implicit minimum-prefix-cost construction,
and inherits neither its contraction/coverage theorem nor its asymptotic
runtime. Pilot work is capped at 4,000,000 block-DP cells per input; an interval
that would exceed the cap uses linear interpolation. All pilot time is included
in reported timings. Equal-length intervals have zero net displacement and
use linear interpolation without probing; this cannot predict compensating
insertion/deletion excursions without informative anchors.

## Exact certificate and proof

Let B contain the source and target of the usual directed edit graph. Compute
g_B(v), the exact shortest prefix cost using only B. Set **U=g_B(target)**,
the restricted optimum; an arbitrary feasible path cost is not sufficient
for the argument below. For each directed edge (u,v) leaving B, use

```
h(v) = abs(remaining_length_first - remaining_length_second)
L = min(U, min_{u in B, v outside B} [g_B(u) + edge_cost(u,v) + h(v)])
```

That is the default `--bound length` heuristic. `--bound composition` and
`--bound qgram` replace h with stronger admissible maxima, defined and proved
in `BOUNDS.md`. `--bound ordered` additionally includes the order-aware
relaxation in `ORDERED_BOUND.md`. `--bound sparse` computes threshold proofs
for that same relaxation as described in `SPARSE_BOUND.md`. The same
certificate argument applies to all five settings.

Only reachable u matter. Insertions, deletions AND diagonal substitution/match
edges are included. An empty exit set has minimum +infinity.

If a globally optimal path stays in B, its cost equals U. Otherwise consider
its first edge leaving B. Its prefix lies in B, so g_B(u) is no larger than
that prefix cost; h(v) is no larger than its remaining cost. That exit's bound
is therefore no larger than the optimal total cost D. In both cases L<=D<=U.

Return exact when L=U. For requested epsilon>=0, returning when
U<=(1+epsilon)*L guarantees a script of cost at most (1+epsilon)*D. The argument
also handles paths leaving and re-entering the band. It does not assume an
anchor is correct or that a particular path touches the boundary.

Otherwise, double the radius around the same proposal and recompute the band.
At radius max(n,m), search the full matrix. Scores use integers throughout;
floating-point prediction can change only which cells are searched, not the
certificate. The default cumulative DP budget is 50,000,000 cells per input;
`--max-cells` changes it. Budget exhaustion is explicitly uncertified. This
budget excludes pilot cells, which have their separate cap, and seed storage.

The engine recomputes expanded bands rather than incrementally updating them
and retains the score matrix for traceback. Even with stronger suffix bounds,
exact certification may require expansion after an optimum has been found.
Worst-case time and storage remain quadratic; some inputs deliberately reach
the full matrix.

## Validation and measurements

`tests/test_correctness.py` uses an independent ordinary full-DP oracle, replays
every returned CIGAR, and checks every finite round's lower/upper interval.
It covers exhaustive small binary strings, random and adversarial sequences,
empty inputs, several starting widths, approximation tolerances, and budget
exhaustion. See `results/correctness.json` for the actual completed run.

`benchmarks/run_bench.py` generates deterministic synthetic pairs, compares
all modes on identical inputs, and records repeated timings and exact-score
agreement. It retains the sequences, raw records, and aggregated summaries.
Defaults use lengths 256, 1024 and 4096, three seeds and three repetitions.
This measures a mechanism in one scalar implementation. It does **not**
establish a speedup over Edlib, A*PA2, QuickEd, WFA2 or BiWFA, and does not
constitute validation on real sequencing data.

`benchmarks/run_bound_bench.py` is the v2 controlled experiment: four unchanged
predictors crossed with all three suffix bounds, the exact original 63 pairs,
55 fresh indel-focused pairs, and five timed repetitions with rotated/reversed
setting order. It records preprocessing and certificate overhead as well as
alignment cells. Forty fresh cases also have an analytically known optimum.
`tests/test_lower_bounds.cpp` independently checks suffix histograms, arbitrary
cursor movement, threshold short-circuiting and admissibility against full DP.
`tests/test_bounds_integration.py` checks complete alignments and certificates,
and compares default/length behavior with the included v1 reference source.

`benchmarks/run_ordered_bench.py` is the v3 experiment: length, qgram and ordered
bounds crossed with the same four predictors, the same 118 pairs, and 30 further
compensated-indel pairs fixed before timing. `tests/test_ordered_bounds.cpp`
checks the optimized index against a slow minimum over every possible target
chunk length and against unrestricted edit distance. The integration audit
compares old modes against the untouched `reference_v2/` source.

`benchmarks/run_sparse_bench.py` is the v4 experiment: the same four predictors
crossed with length, qgram, ordered and sparse bounds, the unchanged 148 pairs,
and 30 new compensated-indel pairs fixed before timing. Five timed repetitions
and one warmup use rotated/reversed setting order. The primary comparison is
total runtime against linear/length, with sparse versus dense reported as a
secondary mechanism check. `tests/test_sparse_bounds.cpp` compares threshold
proofs and cached lower/upper endpoints with a slow all-partition oracle;
`tests/test_sparse_integration.py` audits complete certificates, budget
fallbacks, and old-mode behavior against `reference_v3/`.

`--order-block` defaults to 8 and accepts 1..8. `--order-max-cells` defaults to
4,000,000 and limits the separate relaxed index, including its terminal row.
`--bound-memory` sets the shared bound-index payload cap in bytes (default
67,108,864). Omitting an index because of a limit is reported and preserves
the existing admissible bounds. The main `--max-cells` DP budget is separate.

`--sparse-max-work` defaults to 4,000,000 cumulative nontrivial proof-node
expansions; revisits consume work too. The sparse directory, packed words,
explicit stack and tiles share the same `--bound-memory` payload cap as the
histograms. `--order-block` also applies to sparse mode; `--order-max-cells`
only caps the dense ordered index. Sparse fallbacks retain already proven
bounds, so hitting a bound budget does not itself make an alignment uncertified.

Metrics:

- `dp_cells`: cumulative allocated/visited band cells across all rounds,
  including boundaries and recomputation.
- `probe_cells`: additional block-DP recurrence cells used by prediction.
- `peak_cells`: largest band score allocation; `4*peak_cells` bytes is only
  score storage, not total process memory.
- `predictor_ms`, `dp_ms`, `certificate_ms`, `total_ms`: internal wall-clock
  time; total includes predictor, DP, certificate and traceback, excludes
  output serialization/process startup. Batch process time is saved separately.
- `trace`: each radius, band size, restricted upper and certified lower bound.
  A null upper bound denotes no within-band source-to-target path yet.
- `bound_preprocess_ms`: lazy suffix-index construction, separate from
  `certificate_ms` and included in `total_ms`.
- `bound_bytes`: retained index payload, including gram arrays and cursor
  histograms; excludes allocator metadata and other process memory. The
  default 64MiB index cap is independent of the main DP cell budget.
- `bound_active`, `bound_fallback`: stages actually constructed and any omitted
  for the index cap. An omitted stage weakens the bound and preserves safety.
- `bound_calls`, `composition_calls`, `q2_calls`, `q4_calls`: query work after
  cheaper stages fail to reach the current certificate threshold.
- `bound_cursor_updates`: total coordinate distance traversed across bound
  stages and exit streams, including backward travel.
- `ordered_cells`, `ordered_bytes`, `ordered_calls`: extra relaxed-index cells,
  retained payload and queries. Index work is not part of `dp_cells`.
- `ordered_preprocess_ms`: a subset of `bound_preprocess_ms`, not an additional
  cost to add to it; both are already included in `total_ms`.
- `ordered_fallback`: why the ordered index was omitted, if a limit or allocation
  failure prevented construction. An empty reason can also mean construction
  was never needed; inspect `bound_active` and `ordered_cells`.
- `sparse_queries`, `sparse_expansions`, `sparse_cache_hits`: threshold queries,
  nontrivial expanded nodes (including revisits), and cached interval decisions.
- `sparse_guard_rejections`: cheap feasible upper bounds below the requested
  threshold. This is normal proof avoidance, not a budget fallback.
- `sparse_cached_states`, `sparse_tiles`, `sparse_bytes`, `sparse_peak_stack`:
  touched cache states, 64-column tiles, total sparse payload, and used stack
  depth. The flat tile directory and full-depth stack are included in bytes.
- `sparse_preprocess_ms`: a subset of `bound_preprocess_ms`. Proof traversal
  and tile allocation are in `certificate_ms`; all are in `total_ms`.
- `sparse_fallback`: why new proof work stopped. Existing proven cache values
  remain usable. `+sparse` in `bound_active` means index payload was allocated;
  allocation-free cheap-guard rejections can have zero `sparse_bytes`.

Compare `cost` with `linear`, and **`anchor-cost` with `anchor-linear`** to
isolate the cost component. Comparing `anchor-cost` only with `linear` confounds
the contribution of anchors. Full DP shares the C++ recurrence; independent
correctness evidence comes from the Python oracle.

## Mathematical provenance and baselines

Accessed 9 October 2026. Source URLs point to the evolving release.

- OpenAI family 121, [predicted bands](https://github.com/openai/math/blob/main/preprints/An-Almost-Linear-Approximation-Scheme-for-Edit-Distance-September-24-2026/build/sections/05-bands.tex)
  and [comparison-path coverage](https://github.com/openai/math/blob/main/preprints/An-Almost-Linear-Approximation-Scheme-for-Edit-Distance-September-24-2026/build/sections/06-comparison.tex).
- The independent first-exit certificate is standard shortest-path reasoning;
  it is not claimed as a new theorem from that release.
- [A*PA2, WABI 2024](https://drops.dagstuhl.de/entities/document/10.4230/LIPIcs.WABI.2024.17).
- [QuickEd, Bioinformatics 2025](https://academic.oup.com/bioinformatics/article/41/3/btaf112/8075120).
- [Edlib](https://github.com/Martinsos/edlib), [WFA2](https://github.com/smarco/WFA2-lib).

The package contains an independent experimental implementation, not code
copied from those libraries. Results and known limitations are recorded in
`RESULTS.md`, `RESULTS_V2.md`, `RESULTS_V3.md`, `RESULTS_V4.md`, and `RESULTS_V5.md`.
