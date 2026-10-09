# Stronger suffix bounds in prototype v2

The band predictors are unchanged from v1. This experiment replaces only the
admissible suffix heuristic in the first-exit certificate. Scoring remains
global unit-cost Levenshtein distance, including on the indel-focused dataset.
These inequalities and the certificate are not claimed as new mathematical
results from the OpenAI release.

## Bound and proof

For suffixes A=x[i:] and B=y[j:], let C_q(A,g) count overlapping occurrences of
the q-gram g. Empty or shorter-than-q strings have no q-grams. Define

```
delta_q(g) = C_q(A,g) - C_q(B,g)
P_q = sum_g max(delta_q(g), 0)
N_q = sum_g max(-delta_q(g), 0)
h_q = ceil(max(P_q,N_q) / q)
```

An edit removes at most q old q-grams and introduces at most q new q-grams.
Consequently, one edit can decrease P_q by at most q and N_q by at most q.
Both quantities must reach zero when A has been transformed into B. Any K-edit
transformation therefore satisfies P_q<=qK and N_q<=qK, proving h_q<=ED(A,B).
This argument covers substitutions, insertions and deletions, including edits
near string endpoints.

The available heuristics are:

| CLI setting | Suffix lower bound |
| --- | --- |
| `--bound length` | Absolute remaining-length difference; v1 baseline and default. |
| `--bound composition` | Maximum of length difference and h_1. |
| `--bound qgram` | Maximum of length difference, h_1, h_2 and h_4. |

For q=1, h_1 is the exact edit distance if symbol order is ignored. The maximum
of admissible lower bounds is admissible. Their sum is generally **not** safe:
the same edit can explain discrepancies in several histograms.

For a small alphabet the implementation uses exact gram identities. For larger
alphabets, q2/q4 may merge gram identities into fixed hash buckets. Within each
bucket, opposite-sign discrepancies can cancel; summing positive or negative
parts after merging cannot exceed the original P_q or N_q. Hash collisions can
only weaken the heuristic and cannot create a false certificate.

## Use inside the existing certificate

For restricted prefix costs g_B and restricted optimal score U, compute

```
L = min(U, min over all reachable exits (u,v) of g_B(u)+cost(u,v)+h(v))
```

The first-exit proof in README.md still gives L<=ED(x,y)<=U. Paths that leave
and later re-enter the band are included in that proof. No anchor is required
to occur on an optimal alignment.

When inspecting one exit, the evaluator needs to know whether its candidate
improves the current minimum L. If the length bound or an intermediate
histogram bound already makes the candidate at least L, the remaining stages
can be skipped. Returning the corresponding threshold is safe: the candidate
cannot change the minimum. This optimization does not discard an exit merely
because it looks unlikely, and it does not turn an upper bound into a lower one.

## Incremental evaluation

Moving a suffix start forward by one position removes at most one q-gram: the
one starting at that position. Maintain the difference histogram, its positive
sum and its negative sum. Updating one bin adjusts both sums in constant time.
Moving a cursor backward reverses the same operation.

The certificate uses separate cursor streams for horizontal exits and for
vertical/diagonal exits on each side of the band. This avoids repeatedly jumping
between its left and right sides. Gram identities are precomputed lazily and
reused. The configured q values are fixed; no parameters are selected from
benchmark outcomes. Implementation metrics expose preprocessing time, retained
bound storage and cursor updates so their cost is visible.

For monotone boundaries, cursor work is linear in the coordinate advances and
number of exits per round, up to the constant number of streams and q values.
Backward movement is supported for correctness on other query orders, but its
work depends on total cursor travel. The method does not make arbitrary
nonmonotone queries constant-time.

## Expected strengths and limitations

Length difference can be zero despite separated insertions and deletions.
Composition and q-gram discrepancies can then provide positive lower bounds.
They still ignore much of sequence order and can remain small when inserted
and deleted material has similar word frequencies. Repeated material can make
them especially weak. A higher q does not necessarily dominate a lower q.

Pure insertions or deletions already have a strong whole-sequence length bound;
their performance benefit cannot be assumed from the compensated-indel case.
The benchmark therefore separates insertion, deletion, compensation and repeat
controls, and compares each predictor with itself under the different bounds.
Any observed benefit must include preprocessing and certificate overhead.

See RESULTS_V2.md for measured outcomes and tests. Finite tests support the
implementation; the inequalities above establish admissibility mathematically.
