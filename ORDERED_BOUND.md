# Ordered block and gap relaxation: prototype v3

This bound retains sequence order and charges displacement between source
blocks and their target chunks. It is combined by maximum with the v2 q-gram
bound under `--bound ordered`. Band prediction is unchanged.

The design draws on established seed/gap reasoning, including the ordered
seed and gap-chaining heuristics in Groot Koerkamp and Ivanov,
[Bioinformatics 2024, DOI 10.1093/bioinformatics/btae032](https://academic.oup.com/bioinformatics/article/40/3/btae032/7587511).
Their gap-chaining heuristic uses the maximum of skipped-seed cost and length
displacement, avoiding double-counting the same edits. Our implementation uses
the explicit dense relaxation below. It is not their complete A*PA algorithm,
does not implement match pruning, and is not claimed as a new result from the
OpenAI mathematics release.

## Relaxed block problem

Partition the first sequence into consecutive blocks of at most k characters,
including a final shorter block if necessary. For one source block S and any
target chunk T, assign cost

```
c(S,T) = 0                         if S == T
         max(1, abs(len(S)-len(T))) otherwise
```

This is a lower bound on their real edit distance: unequal strings require at
least one edit, and length imbalance requires at least that many indels.
The maximum combines those obligations without charging an edit twice.

Let source block boundaries be a_0=0,...,a_B=n. At each boundary define

```
F_B(j) = m-j
F_b(j) = min over t in [j,m] of
         c(x[a_b:a_(b+1)], y[j:t]) + F_(b+1)(t)
```

Target chunks occur in order and may be empty. The terminal row charges target
characters left over after the source is consumed. Every possible exact seed
match is represented implicitly; no unique-match assumption, selected anchor
chain or discarded repetitive hit is required.

To prove admissibility, take any real alignment and cut it at the first crossing
of each source block boundary. Its pieces have disjoint edges; an insertion at
a boundary belongs to one piece only. Each piece's real cost is at least its
relaxed cost, and the terminal insertion cost is charged separately. Minimizing
over all ordered target partitions therefore cannot exceed the real optimum.

## Three-candidate recurrence

Each F row is 1-Lipschitz in target position:

```
abs(F_b(j+1)-F_b(j)) <= 1
```

Adding or removing one target character changes the relaxed cost of one chunk
by at most one. Apply that operation to an optimal partition (or its terminal
insertion tail) in each direction to obtain the inequality. This includes the
special exact-match case: insertion/deletion changes its cost from zero to one.

Let G=F_(b+1), l=a_(b+1)-a_b and p=j+l. For p<=m, the full minimum is

```
F_b(j) = min(G(p) + [source_block != y[j:p]],
             1 + G(p-1),
             1 + G(p+1))
```

Omit the last candidate if p=m. The left neighbor is valid because l>=1.
The Lipschitz property guarantees that among all shorter target chunks the
nearest shorter length attains a minimum, and likewise for longer chunks.
When p>m, the recurrence is simply

```
F_b(j) = p-m + G(m)
```

Target words of lengths up to eight bytes have collision-free packed 64-bit
identities. The implementation precomputes identities for the full block size
and, when needed, the final shorter block size. Equality checks therefore take
constant time even in repeats. `--order-block` accepts 1 through 8 and defaults
to 8; the benchmark fixes it at 8 before inspecting results.

## Queries inside a source block

A certificate exit can start anywhere, not only at a block boundary. Let s be
the next source boundary at or after i, and r=s-i. Ignore the sequence content
of this partial source block, retaining only its length obligation. Then

```
h_order(i,j) = F_s(min(m,j+r)) + max(0,j+r-m)
```

Here F_s denotes the row for boundary s. This formula equals the minimum over
all target-prefix lengths t>=0 of `abs(r-t)+F_s(j+t)`, by the same Lipschitz
property. It is therefore a lower bound, rather than a choice of one supposedly
correct target crossing. At a block boundary r=0 and the exact relaxed row
value is used. Do not add an unmatched-seed charge for the ignored partial block.

## Integration and cost

`--bound ordered` first checks the cheaper length/composition/q-gram bounds.
Only if they cannot rule out improving the current certificate minimum does it
build or query this index. The returned suffix bound is a maximum, never a sum.
The independent first-exit proof in README.md is unchanged.

There are `(ceil(n/k)+1)*(m+1)` index cells, including the terminal row.
Preprocessing and storage are O(nm/k), with constant-time suffix queries after
construction. For fixed k this remains quadratic preprocessing; it is a trade
of index work for fewer expanded alignment bands, not an asymptotic breakthrough.

The index is constructed at most once per input and reused across rounds.
`--order-max-cells` defaults to 4,000,000. The existing 64MiB bound-payload limit
also applies. Exceeding a limit omits this component and records a fallback to
the v2 bounds; it must never certify using a truncated or incomplete index.
Index cells, bytes, preprocessing time, queries and fallback reasons are reported.
Preprocessing is included in total runtime and must be counted when comparing
work with earlier variants.

## Limits of the relaxation

Different equal-length blocks can cost only one, even if every character differs.
It is not exact edit distance even at k=1: `A` versus `CC` has relaxed cost 1
but real edit distance 2. This is deliberate underestimation, not a correctness
defect.

At the starting state of equal-length strings, pairing same-length chunks gives
an upper limit of ceil(n/k) on the relaxed bound. For n=4096 and k=8 it therefore
cannot exceed 512, even if real distance is 818. Sequence order helps, but this
ceiling can still prevent a narrow-band exact certificate. The v2 histogram
components can sometimes be stronger, which is why the maximum is retained.

The independent tests compare the optimized recurrence against a slow minimum
over every target chunk length, and compare arbitrary suffix queries against
unrestricted edit distance. RESULTS_V3.md records actual tests and measurements.
