# Sparse threshold certification: prototype v4

The sparse evaluator targets the **same relaxed ordered-block value F** as v3.
The band predictors and scoring model are unchanged. This experiment changes
how much of the bound is evaluated and cached; it does not introduce a stronger
mathematical relaxation. See `ORDERED_BOUND.md` for the definition and proof of F.

`--bound sparse` combines the v2 length/composition/q-gram bounds with this
evaluator. All components are combined by maximum. The independent first-exit
certificate in README.md remains responsible for exactness.

## The threshold question

At a band exit, suppose the prefix and edge have cost p and the current minimum
certificate value is L. This exit cannot lower that minimum if its suffix bound
reaches T=L-p. There is no need to compute the full relaxed value in that case.

The v3 recurrence has at most three candidates:

```
F(state) = min over children r of cost_r + F(child_r)
```

Consequently, proving `F(state) >= T` requires proving **every** child has
`F(child_r) >= T-cost_r`. One unsuccessful child is enough to stop that proof
attempt. The sparse implementation follows these dependencies with an explicit
stack and reuses bounds from earlier queries and expansion rounds.

## Cache intervals, not guesses

Each cached state has a proven lower and upper bound on its true relaxed F.
Initially, the remaining-length difference is a valid lower bound. The upper
bound below supplies additional inexpensive information.

- A lower bound reaching T proves success immediately.
- An upper bound below T proves the threshold is unreachable.
- If every child proves its required threshold, the parent's lower bound can
  be raised to T.
- If a child proves its value is below its required threshold, the parent's
  relaxed value is strictly below T. Integer costs allow its upper bound to be
  lowered to T-1.
- Work or memory exhaustion is **unknown**, not a disproved threshold. It may
  neither raise a lower bound nor lower an upper bound.

The value returned to the alignment certificate is always a proven lower bound.
Upper bounds are used only to avoid futile proof attempts. A minimum over only
the children visited so far is not a valid lower bound and is never returned.

## A cheap upper bound that avoids futile traversal

At a full source-block boundary, let R be the number of remaining blocks and
let delta be the difference between remaining source and target lengths. Then

```
F <= max(R, abs(delta))
```

To construct a relaxed partition achieving this upper bound, distribute the
net length discrepancy across the R chunks. When its magnitude is at most R,
adjust each target chunk's length by at most one; each chunk costs at
most one, including unchanged-length chunks. When its magnitude exceeds R,
assign at least one unit of length change to every chunk and distribute the
rest; the summed relaxed costs equal the discrepancy. For deletions, capacity
is sufficient because every block has positive length and their total source
length covers the required deletion. With R=0, only terminal insertions remain.

If this upper bound is below T, even the full dense ordered index could not
close that exit at the requested threshold. The traversal can be skipped
without allocating the score tiles. This is a property of the relaxation,
not a statistical prediction about the real alignment.

## Partial source blocks

For a query at source position i, let s be the next source-block boundary,
r=s-i, j'=min(m,j+r), and o=max(0,j+r-m). The v3 relaxed query is

```
h_order(i,j) = F_s(j') + o
```

The sparse evaluator asks whether `F_s(j') >= T-o`. On success it may return T;
otherwise it returns a valid cached or basic lower bound at that state plus o.
The true relaxed rows are 1-Lipschitz, as proved for v3. Partially populated cache
rows need not be Lipschitz; no such property is assumed about them.

## Sparse storage and bounded work

Scores are cached in 64-column tiles, allocated only when needed. A direct
pointer directory supports constant-time tile lookup. The directory itself has
O(ceil(n/k)*ceil((m+1)/64)) entries: this is not a fully sparse asymptotic data
structure, and its memory is checked before allocation. Exact packed word
identities are prepared lazily in O(n+m) work. The explicit proof stack avoids
dependence on the process recursion limit.

`--sparse-max-work` defaults to 4,000,000 nontrivial proof-node expansions per
input. Revisits count against this budget; caching does not imply each node is
expanded only once. `--bound-memory` is the shared payload cap for histogram
indices, the sparse directory, word identities, stack and allocated tiles.
Allocator bookkeeping and other process memory are excluded from this payload
metric. Limit exhaustion is reported and leaves the existing admissible bounds
available. No incomplete proof is treated as a success.

Metrics expose proof queries, expanded nodes, cache hits, cheap-guard rejections,
cached states, tiles, payload bytes, stack use and preprocessing time. Expanded
proof nodes and dense DP cells are different operations; total runtime is the
primary performance comparison. Sparse preprocessing is part of total bound
preprocessing, and proof traversal is included in certificate time.

## Relationship to the dense version

With sufficient budgets, the proof procedure answers `F>=T` exactly. In exact
alignment mode it therefore preserves the dense ordered certificate's stopping
round: if every exit can reach U, every corresponding proof succeeds; if one
cannot, dense certification also fails. Intermediate lower bounds on failed
rounds can be weaker, because the sparse evaluator need not compute the full F.

With budget fallbacks, the sparse component can be weaker than dense F but
remains admissible. Positive-epsilon runs retain the approximation guarantee;
they need not return at the same round or with the same approximate score as
dense evaluation. The independent integration tests check these distinctions.

Worst-case work can still approach dense evaluation, and threshold revisits add
overhead. Sparsity alone is not evidence of a speedup. RESULTS_V4.md reports the
fixed experiment against ordinary linear banding as well as earlier bounds.
